"use strict";
/* global Ammo,THREE */

const Type = (exports.Type = {
  BOX: "box",
  CYLINDER: "cylinder",
  SPHERE: "sphere",
  CAPSULE: "capsule",
  CONE: "cone",
  HULL: "hull",
  MESH: "mesh"
});

exports.createCollisionShape = (function() {
  const halfExtents = new THREE.Vector3();
  const pos = new THREE.Vector3();
  const quat = new THREE.Quaternion();
  const scale = new THREE.Vector3();
  const box = new THREE.Box3();

  return function(sceneRoot, options) {
    const autoGenerateShape = options.hasOwnProperty("autoGenerateShape") ? options.autoGenerateShape : true;
    const mergeGeometry = options.hasOwnProperty("mergeGeometry") ? options.mergeGeometry : true;
    const type = options.type || Type.HULL;
    const recenter = options.hasOwnProperty("recenter") ? options.recenter : false;
    if (options.halfExtents) {
      halfExtents.set(options.halfExtents.x, options.halfExtents.y, options.halfExtents.z);
    }
    const minHalfExtent = options.hasOwnProperty("minHalfExtent") ? options.minHalfExtent : 0;
    const maxHalfExtent = options.hasOwnProperty("maxHalfExtent") ? options.maxHalfExtent : Number.POSITIVE_INFINITY;
    const cylinderAxis = options.cylinderAxis || "y";
    const sphereRadius = options.sphereRadius || NaN;
    const margin = options.hasOwnProperty("margin") ? options.margin : 0.01;
    const hullMaxVertices = options.hullMaxVertices || 100000;

    let collisionShape;
    let triMesh;
    let shapeHull;

    let meshes;
    let vertices;

    if ((mergeGeometry || autoGenerateShape || recenter) && !sceneRoot) {
      console.warn("cannot use mergeGeometry, autoGenerateShape, or recenter if sceneRoot is null");
      return;
    }

    if (mergeGeometry) {
      meshes = _getMeshes(sceneRoot);
    } else {
      meshes = [sceneRoot];
    }

    if (type !== "mesh") {
      if (recenter) {
        _recenter(sceneRoot, meshes);
      }
    }

    if (autoGenerateShape) {
      vertices = _getVertices(sceneRoot, meshes);

      if ([Type.SPHERE, Type.HULL, Type.MESH].indexOf(type) === -1) {
        box.setFromPoints(vertices);
        const { max, min } = box;
        halfExtents
          .subVectors(max, min)
          .multiplyScalar(0.5)
          .clampScalar(minHalfExtent, maxHalfExtent);
      }
    }
    const { x, y, z } = halfExtents;

    //TODO: Support convex hull decomposition, compound shapes, gimpact (dynamic trimesh)
    switch (type) {
      case Type.BOX: {
        const halfExtents = new Ammo.btVector3(x, y, z);
        collisionShape = new Ammo.btBoxShape(halfExtents);
        Ammo.destroy(halfExtents);
        break;
      }
      case Type.SPHERE: {
        let radius = 1;
        if (sphereRadius) {
          radius = sphereRadius;
        } else {
          const sphere = new THREE.Sphere();
          sphere.setFromPoints(vertices);
          if (isFinite(sphere.radius)) {
            radius = sphere.radius;
          }
        }
        collisionShape = new Ammo.btSphereShape(radius);
        break;
      }
      case Type.CYLINDER: {
        const halfExtents = new Ammo.btVector3(x, y, z);
        switch (cylinderAxis) {
          case "y":
            collisionShape = new Ammo.btCylinderShape(halfExtents);
            break;
          case "x":
            collisionShape = new Ammo.btCylinderShapeX(halfExtents);
            break;
          case "z":
            collisionShape = new Ammo.btCylinderShapeZ(halfExtents);
            break;
        }
        Ammo.destroy(halfExtents);
        break;
      }
      case Type.CAPSULE: {
        switch (cylinderAxis) {
          case "y":
            collisionShape = new Ammo.btCapsuleShape(Math.max(x, z), y * 2);
            break;
          case "x":
            collisionShape = new Ammo.btCapsuleShapeX(Math.max(y, z), x * 2);
            break;
          case "z":
            collisionShape = new Ammo.btCapsuleShapeZ(Math.max(x, y), z * 2);
            break;
        }
        break;
      }
      case Type.CONE: {
        switch (cylinderAxis) {
          case "y":
            collisionShape = new Ammo.btConeShape(Math.max(x, z), y * 2);
            break;
          case "x":
            collisionShape = new Ammo.btConeShapeX(Math.max(y, z), x * 2);
            break;
          case "z":
            collisionShape = new Ammo.btConeShapeZ(Math.max(x, y), z * 2);
            break;
        }
        break;
      }
      case Type.HULL: {
        sceneRoot.matrixWorld.decompose(pos, quat, scale);
        const localScale = new Ammo.btVector3(scale.x, scale.y, scale.z);
        const vec3 = new Ammo.btVector3();
        const originalHull = new Ammo.btConvexHullShape();
        originalHull.setMargin(margin);
        if (vertices.length > hullMaxVertices) {
          console.warn(
            "too many vertices for hull shape; randomly sampling " +
              hullMaxVertices +
              " from " +
              vertices.length +
              " vertices"
          );
          vertices = getRandomSample(vertices, hullMaxVertices);
        }
        for (let i = 0; i < vertices.length; i++) {
          vec3.setValue(vertices[i].x, vertices[i].y, vertices[i].z);
          originalHull.addPoint(vec3, i == vertices.length - 1);
        }

        collisionShape = originalHull;
        if (originalHull.getNumVertices() >= 100 || true) {
          //Bullet documentation says don't use convexHulls with 100 verts or more
          shapeHull = new Ammo.btShapeHull(originalHull);
          shapeHull.buildHull(margin);
          Ammo.destroy(originalHull);
          collisionShape = new Ammo.btConvexHullShape(
            Ammo.getPointer(shapeHull.getVertexPointer()),
            shapeHull.numVertices()
          );
        }
        collisionShape.setLocalScaling(localScale);

        Ammo.destroy(localScale);
        Ammo.destroy(vec3);
        break;
      }
      case Type.MESH: {
        sceneRoot.matrixWorld.decompose(pos, quat, scale);
        const localScale = new Ammo.btVector3(scale.x, scale.y, scale.z);
        const a = new Ammo.btVector3();
        const b = new Ammo.btVector3();
        const c = new Ammo.btVector3();
        triMesh = new Ammo.btTriangleMesh(true, false);

        for (let j = 0; j < vertices.length; j += 3) {
          a.setValue(vertices[j].x, vertices[j].y, vertices[j].z);
          b.setValue(vertices[j + 1].x, vertices[j + 1].y, vertices[j + 1].z);
          c.setValue(vertices[j + 2].x, vertices[j + 2].y, vertices[j + 2].z);
          triMesh.addTriangle(a, b, c, j == vertices.length - 3);
          //TODO: limit number of triangles?
        }

        collisionShape = new Ammo.btBvhTriangleMeshShape(triMesh, true, true);
        collisionShape.setLocalScaling(localScale);

        Ammo.destroy(localScale);
        Ammo.destroy(a);
        Ammo.destroy(b);
        Ammo.destroy(c);
        break;
      }

      default:
        console.warn(shape + " is not currently supported");
        return;
    }

    collisionShape.setMargin(margin);

    collisionShape.destroy = () => {
      if (shapeHull) Ammo.destroy(shapeHull);
      if (triMesh) Ammo.destroy(triMesh);
      if (collisionShape) Ammo.destroy(collisionShape);
    };

    collisionShape.type = type;

    return collisionShape;
  };
})();

//https://stackoverflow.com/a/37835673
const getRandomSample = (function() {
  const swaps = [];
  return function(array, size) {
    let r,
      i = array.length,
      end = i - size,
      temp;

    while (i-- > end) {
      r = Math.floor(Math.random() * (i + 1));
      temp = array[r];
      array[r] = array[i];
      array[i] = temp;
      swaps.push(i);
      swaps.push(r);
    }

    const sample = array.slice(end);

    while (size--) {
      i = swaps.pop();
      r = swaps.pop();
      temp = array[i];
      array[i] = array[r];
      array[r] = temp;
    }

    return sample;
  };
})();


// Whether this THREE.Object3D is a mesh we should take into account for our shape.
function _shouldInclude(obj) {
  return obj.isMesh && (!THREE.Sky || obj.__proto__ != THREE.Sky.prototype);
}

function _getMeshes(sceneRoot) {
  let meshes = [];
  sceneRoot.traverse(o => {
    if (_shouldInclude(o)) {
      meshes.push(o);
    }
  });
  return meshes;
}

const _getVertices = (function() {
  const vertexPool = [];
  const vertices = [];
  const matrix = new THREE.Matrix4();
  const inverse = new THREE.Matrix4();

  return function(sceneRoot, meshes) {
    while (vertices.length > 0) {
      vertexPool.push(vertices.pop());
    }

    inverse.getInverse(sceneRoot.matrixWorld);

    for (let j = 0; j < meshes.length; j++) {
      const mesh = meshes[j];

      const geometry = mesh.geometry.index ? mesh.geometry.toNonIndexed() : mesh.geometry.clone();

      if (mesh !== sceneRoot) {
        matrix.multiplyMatrices(inverse, mesh.matrixWorld);
        geometry.applyMatrix(matrix);
      }

      if (geometry.isBufferGeometry) {
        const components = geometry.attributes.position.array;
        for (let i = 0; i < components.length; i += 3) {
          const x = components[i];
          const y = components[i + 1];
          const z = components[i + 2];

          if (vertexPool.length > 0) {
            vertices.push(vertexPool.pop().set(x, y, z));
          } else {
            vertices.push(new THREE.Vector3(x, y, z));
          }
        }
      } else {
        for (let i = 0; i < geometry.vertices.length; i++) {
          const vertex = geometry.vertices[i];
          if (vertexPool.length > 0) {
            vertices.push(vertexPool.pop().copy(vertex));
          } else {
            vertices.push(new THREE.Vector3(vertex.x, vertex.y, vertex.z));
          }
        }
      }
    }

    return vertices;
  };
})();

const _recenter = (function() {
  const geometries = [];
  const offset = new THREE.Matrix4();
  const center = new THREE.Vector3();

  return function(sceneRoot, meshes) {
    if (meshes.length === 1) {
      meshes[0].geometry.center();
      return;
    }

    const { min, max } = _getBoundingBox(meshes);
    center.addVectors(max, min).multiplyScalar(-0.5);
    offset.makeTranslation(center.x, center.y, center.z);

    for (let j = 0; j < meshes.length; j++) {
      const mesh = meshes[j];
      if (geometries.indexOf(mesh.geometry.uuid) !== -1) {
        continue;
      }
      mesh.geometry.applyMatrix(offset);
      geometries.push(mesh.geometry.uuid);
    }
  };
})();

const _getBoundingBox = (function() {
  const boundingBox = {
    min: new THREE.Vector3(Number.MAX_VALUE),
    max: new THREE.Vector3(Number.MIN_VALUE)
  };

  return function(meshes) {
    for (let i = 0; i < meshes.length; ++i) {
      const mesh = meshes[i];
      if (!mesh.geometry.boundingBox) {
        mesh.geometry.computeBoundingBox();
      }
      const box = mesh.geometry.boundingBox;

      boundingBox.min.x = Math.min(box.min.x, box.min.x);
      boundingBox.min.y = Math.min(box.min.y, box.min.y);
      boundingBox.min.z = Math.min(box.min.z, box.min.z);

      boundingBox.max.x = Math.max(box.max.x, box.max.x);
      boundingBox.max.y = Math.max(box.max.y, box.max.y);
      boundingBox.max.z = Math.max(box.max.z, box.max.z);
    }
    return boundingBox;
  };
})();
