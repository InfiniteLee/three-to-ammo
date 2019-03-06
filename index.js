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

  return function(sceneRoot, options) {
    const autoGenerateShape = options.hasOwnProperty("autoGenerateShape") ? options.autoGenerateShape : true;
    const mergeGeometry = options.hasOwnProperty("mergeGeometry") ? options.mergeGeometry : true;
    const type = options.type || Type.HULL;
    const recenter = options.hasOwnProperty("recenter") ? options.recenter : false;
    const minHalfExtent = options.hasOwnProperty("minHalfExtent") ? options.minHalfExtent : 0;
    const maxHalfExtent = options.hasOwnProperty("maxHalfExtent") ? options.maxHalfExtent : Number.POSITIVE_INFINITY;
    const cylinderAxis = options.cylinderAxis || "y";
    const margin = options.hasOwnProperty("margin") ? options.margin : 0.01;
    const hullMaxVertices = options.hullMaxVertices || 100000;

    let sphereRadius = options.sphereRadius || NaN;
    let collisionShape;
    let triMesh;
    let shapeHull;

    let meshes;
    let vertices;

    if ((mergeGeometry || autoGenerateShape || recenter) && !sceneRoot) {
      console.warn("cannot use mergeGeometry, autoGenerateShape, or recenter if sceneRoot is null");
      return null;
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

    //TODO: Support convex hull decomposition, compound shapes, gimpact (dynamic trimesh)
    switch (type) {
    case Type.BOX: {
      const bounds = options.halfExtents || computeHalfExtents(sceneRoot, meshes, halfExtents).clampScalar(minHalfExtent, maxHalfExtent);
      collisionShape = _createBoxShape(bounds);
      break;
    }
    case Type.CYLINDER: {
      const bounds = options.halfExtents || computeHalfExtents(sceneRoot, meshes, halfExtents).clampScalar(minHalfExtent, maxHalfExtent);
      collisionShape = _createCylinderShape(bounds, cylinderAxis);
      break;
    }
    case Type.CAPSULE: {
      const bounds = options.halfExtents || computeHalfExtents(sceneRoot, meshes, halfExtents).clampScalar(minHalfExtent, maxHalfExtent);
      collisionShape = _createCapsuleShape(bounds, cylinderAxis);
      break;
    }
    case Type.CONE: {
      const bounds = options.halfExtents || computeHalfExtents(sceneRoot, meshes, halfExtents).clampScalar(minHalfExtent, maxHalfExtent);
      collisionShape = _createConeShape(bounds, cylinderAxis);
      break;
    }
    case Type.SPHERE: {
      const radius = options.sphereRadius || computeRadius(_getVertices(sceneRoot, meshes));
      collisionShape = new Ammo.btSphereShape(radius);
      break;
    }
    case Type.HULL: {
      let vertices = _getVertices(sceneRoot, meshes);
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
      collisionShape = _createHullShape(sceneRoot, vertices, margin);
      break;
    }
    case Type.MESH: {
      collisionShape = _createTriMeshShape(sceneRoot, _getVertices(sceneRoot, meshes));
      break;
    }
    default:
      console.warn(type + " is not currently supported");
      return null;
    }

    collisionShape.type = type;
    collisionShape.setMargin(margin);
    collisionShape.destroy = () => {
      for (let res of collisionShape.resources || []) {
        Ammo.destroy(res);
      }
      Ammo.destroy(collisionShape);
    };

    return collisionShape;
  };
})();

const _createBoxShape = function({x, y, z} = halfExtents) {
  const btHalfExtents = new Ammo.btVector3(x, y, z);
  const collisionShape = new Ammo.btBoxShape(btHalfExtents);
  Ammo.destroy(btHalfExtents);
  return collisionShape;
};

const _createCylinderShape = function({x, y, z} = halfExtents, cylinderAxis) {
  const btHalfExtents = new Ammo.btVector3(x, y, z);
  const collisionShape = (() => {
    switch (cylinderAxis) {
    case "y":
      return new Ammo.btCylinderShape(btHalfExtents);
    case "x":
      return new Ammo.btCylinderShapeX(btHalfExtents);
    case "z":
      return new Ammo.btCylinderShapeZ(btHalfExtents);
    }
    return null;
  })();
  Ammo.destroy(btHalfExtents);
  return collisionShape;
};

const _createConeShape = function({x, y, z} = halfExtents, cylinderAxis) {
  switch (cylinderAxis) {
  case "y":
    return new Ammo.btConeShape(Math.max(x, z), y * 2);
  case "x":
    return new Ammo.btConeShapeX(Math.max(y, z), x * 2);
  case "z":
    return new Ammo.btConeShapeZ(Math.max(x, y), z * 2);
  }
  return null;
};

const _createCapsuleShape = function({x, y, z} = halfExtents, capsuleAxis) {
  switch (capsuleAxis) {
  case "y":
    return new Ammo.btCapsuleShape(Math.max(x, z), y * 2);
  case "x":
    return new Ammo.btCapsuleShapeX(Math.max(y, z), x * 2);
  case "z":
    return new Ammo.btCapsuleShapeZ(Math.max(x, y), z * 2);
  }
  return null;
};

const _createHullShape = (function() {
  const pos = new THREE.Vector3();
  const quat = new THREE.Quaternion();
  const scale = new THREE.Vector3();
  return function(sceneRoot, vertices, margin) {
    sceneRoot.matrixWorld.decompose(pos, quat, scale);
    const localScale = new Ammo.btVector3(scale.x, scale.y, scale.z);
    const vec3 = new Ammo.btVector3();
    const originalHull = new Ammo.btConvexHullShape();
    originalHull.setMargin(margin);
    for (let i = 0; i < vertices.length; i++) {
      vec3.setValue(vertices[i].x, vertices[i].y, vertices[i].z);
      originalHull.addPoint(vec3, i == vertices.length - 1);
    }
    let collisionShape = originalHull;
    if (originalHull.getNumVertices() >= 100 || true) {
      //Bullet documentation says don't use convexHulls with 100 verts or more
      const shapeHull = new Ammo.btShapeHull(originalHull);
      shapeHull.buildHull(margin);
      Ammo.destroy(originalHull);
      collisionShape = new Ammo.btConvexHullShape(
        Ammo.getPointer(shapeHull.getVertexPointer()),
        shapeHull.numVertices()
      );
      collisionShape.resources = [shapeHull];
    }
    collisionShape.setLocalScaling(localScale);
    Ammo.destroy(localScale);
    Ammo.destroy(vec3);
    return collisionShape;
  };
})();

const _createTriMeshShape = (function() {
  const pos = new THREE.Vector3();
  const quat = new THREE.Quaternion();
  const scale = new THREE.Vector3();
  return function(sceneRoot, vertices) {
    sceneRoot.matrixWorld.decompose(pos, quat, scale);
    const localScale = new Ammo.btVector3(scale.x, scale.y, scale.z);
    const a = new Ammo.btVector3();
    const b = new Ammo.btVector3();
    const c = new Ammo.btVector3();
    const triMesh = new Ammo.btTriangleMesh(true, false);

    for (let j = 0; j < vertices.length; j += 3) {
      a.setValue(vertices[j].x, vertices[j].y, vertices[j].z);
      b.setValue(vertices[j + 1].x, vertices[j + 1].y, vertices[j + 1].z);
      c.setValue(vertices[j + 2].x, vertices[j + 2].y, vertices[j + 2].z);
      triMesh.addTriangle(a, b, c, j == vertices.length - 3);
      //TODO: limit number of triangles?
    }

    const collisionShape = new Ammo.btBvhTriangleMeshShape(triMesh, true, true);
    collisionShape.setLocalScaling(localScale);
    collisionShape.resources = [triMesh];
    Ammo.destroy(localScale);
    Ammo.destroy(a);
    Ammo.destroy(b);
    Ammo.destroy(c);
    return collisionShape;
  };
})();

const computeRadius = (function() {
  const sphere = new THREE.Sphere();
  return function(vertices) {
    sphere.setFromPoints(vertices);
    return isFinite(sphere.radius) ? sphere.radius : 1;
  };
})();

const getBufferGeometry = (function() {
  const bufferGeometry = new THREE.BufferGeometry();
  return function(geo) {
    return geo.isBufferGeometry ? geo : bufferGeometry.fromGeometry(geo);
  };
})();

const computeHalfExtents = (function() {
  const v = new THREE.Vector3();
  const matrix = new THREE.Matrix4();
  const inverse = new THREE.Matrix4();
  const min = new THREE.Vector3();
  const max = new THREE.Vector3();
  const bufferGeometry = new THREE.BufferGeometry();
  return function(sceneRoot, meshes, target) {

    let minX = + Infinity;
    let minY = + Infinity;
    let minZ = + Infinity;
    let maxX = - Infinity;
    let maxY = - Infinity;
    let maxZ = - Infinity;

    inverse.getInverse(sceneRoot.matrixWorld);

    for (let mesh of meshes) {
      const geometry = getBufferGeometry(mesh.geometry);
      const components = geometry.attributes.position.array;
      if (mesh !== sceneRoot) {
        matrix.multiplyMatrices(inverse, mesh.matrixWorld);
      } else {
        matrix.identity();
      }
      for (let i = 0; i < components.length; i += 3) {
        v.set(components[i], components[i + 1], components[i + 2]);
        v.applyMatrix4(matrix);
		if ( v.x < minX ) minX = v.x;
		if ( v.y < minY ) minY = v.y;
		if ( v.z < minZ ) minZ = v.z;
		if ( v.x > maxX ) maxX = v.x;
		if ( v.y > maxY ) maxY = v.y;
		if ( v.z > maxZ ) maxZ = v.z;
      }
    }

    return target.set(maxX - minX, maxY - minY, maxZ - minZ).multiplyScalar(0.5);
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
