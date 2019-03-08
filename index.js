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
  const bounds = new THREE.Box3();
  const sphere = new THREE.Sphere();
  const halfExtents = new THREE.Vector3();

  return function(sceneRoot, options) {
    const autoGenerateShape = options.hasOwnProperty("autoGenerateShape") ? options.autoGenerateShape : true;
    const mergeGeometry = options.hasOwnProperty("mergeGeometry") ? options.mergeGeometry : true;
    const type = options.type || Type.HULL;
    const minHalfExtent = options.hasOwnProperty("minHalfExtent") ? options.minHalfExtent : 0;
    const maxHalfExtent = options.hasOwnProperty("maxHalfExtent") ? options.maxHalfExtent : Number.POSITIVE_INFINITY;
    const cylinderAxis = options.cylinderAxis || "y";
    const margin = options.hasOwnProperty("margin") ? options.margin : 0.01;

    if ((mergeGeometry || autoGenerateShape) && !sceneRoot) {
      console.warn("cannot use mergeGeometry or autoGenerateShape if sceneRoot is null");
      return null;
    }

    const computeRadius = function() {
      _computeBounds(sceneRoot, mergeGeometry, bounds);
      _computeSphere(sceneRoot, mergeGeometry, bounds, sphere);
      return sphere.radius;
    };

    const computeHalfExtents = function() {
      _computeBounds(sceneRoot, mergeGeometry, bounds);
      return halfExtents
        .subVectors(bounds.max, bounds.min)
        .multiplyScalar(0.5)
        .clampScalar(minHalfExtent, maxHalfExtent);
    };

    //TODO: Support convex hull decomposition, compound shapes, gimpact (dynamic trimesh)
    let collisionShape;
    switch (type) {
      case Type.BOX: {
        const hx = autoGenerateShape ? computeHalfExtents() : options.halfExtents;
        collisionShape = _createBoxShape(hx);
        break;
      }
      case Type.CYLINDER: {
        const hx = autoGenerateShape ? computeHalfExtents() : options.halfExtents;
        collisionShape = _createCylinderShape(hx, cylinderAxis);
        break;
      }
      case Type.CAPSULE: {
        const hx = autoGenerateShape ? computeHalfExtents() : options.halfExtents;
        collisionShape = _createCapsuleShape(hx, cylinderAxis);
        break;
      }
      case Type.CONE: {
        const hx = autoGenerateShape ? computeHalfExtents() : options.halfExtents;
        collisionShape = _createConeShape(hx, cylinderAxis);
        break;
      }
      case Type.SPHERE: {
        const radius = !isNaN(options.sphereRadius) ? options.sphereRadius : computeRadius();
        collisionShape = new Ammo.btSphereShape(radius);
        break;
      }
      case Type.HULL: {
        collisionShape = _createHullShape(sceneRoot, mergeGeometry, margin, options.hullMaxVertices || 100000);
        break;
      }
      case Type.MESH: {
        collisionShape = _createTriMeshShape(sceneRoot, mergeGeometry);
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

// Calls `cb(mesh)` for each mesh under `root` whose vertices we should take into account for the physics shape.
function _iterateMeshes(root, mergeGeometry, cb) {
  if (!mergeGeometry) {
    cb(root);
  } else {
    root.traverse(obj => {
      if (obj.isMesh && (!THREE.Sky || obj.__proto__ != THREE.Sky.prototype)) {
        cb(obj);
      }
    });
  }
}

// Calls `cb(geo, xform)` for each geometry under `root` whose vertices we should take into account for the physics shape.
// `xform` is the transform required to transform the given geometry's vertices into root-local space.
const _iterateGeometries = (function() {
  const transform = new THREE.Matrix4();
  const inverse = new THREE.Matrix4();
  const bufferGeometry = new THREE.BufferGeometry();
  return function(root, mergeGeometry, cb) {
    inverse.getInverse(root.matrixWorld);
    _iterateMeshes(root, mergeGeometry, mesh => {
      const geo = mesh.geometry.isBufferGeometry ? mesh.geometry : bufferGeometry.fromGeometry(mesh.geometry);
      if (mesh !== root) {
        transform.multiplyMatrices(inverse, mesh.matrixWorld);
      } else {
        transform.identity();
      }
      // todo: might want to return null xform if this is the root so that callers can avoid multiplying
      // things by the identity matrix
      cb(geo, transform);
    });
  };
})();

// Sets `target` to the bounding sphere for the geometries underneath `root`.
const _computeSphere = (function() {
  const v = new THREE.Vector3();
  return function(root, mergeGeometry, bounds, target) {
    let maxRadiusSq = 0;
    let { x: cx, y: cy, z: cz } = bounds.getCenter(target.center);
    _iterateGeometries(root, mergeGeometry, (geo, transform) => {
      const components = geo.attributes.position.array;
      for (let i = 0; i < components.length; i += 3) {
        v.set(components[i], components[i + 1], components[i + 2]).applyMatrix4(transform);
        const dx = cx - v.x;
        const dy = cy - v.y;
        const dz = cz - v.z;
        maxRadiusSq = Math.max(maxRadiusSq, dx * dx + dy * dy + dz * dz);
      }
    });
    target.radius = Math.sqrt(maxRadiusSq);
    return target;
  };
})();

// Sets `target` to the bounding box for the geometries underneath `root`.
const _computeBounds = (function() {
  const v = new THREE.Vector3();
  return function(root, mergeGeometry, target) {
    let minX = +Infinity;
    let minY = +Infinity;
    let minZ = +Infinity;
    let maxX = -Infinity;
    let maxY = -Infinity;
    let maxZ = -Infinity;
    _iterateGeometries(root, mergeGeometry, (geo, transform) => {
      const components = geo.attributes.position.array;
      for (let i = 0; i < components.length; i += 3) {
        v.set(components[i], components[i + 1], components[i + 2]).applyMatrix4(transform);
        if (v.x < minX) minX = v.x;
        if (v.y < minY) minY = v.y;
        if (v.z < minZ) minZ = v.z;
        if (v.x > maxX) maxX = v.x;
        if (v.y > maxY) maxY = v.y;
        if (v.z > maxZ) maxZ = v.z;
      }
    });
    target.min.set(minX, minY, minZ);
    target.max.set(maxX, maxY, maxZ);
    return target;
  };
})();

const _createBoxShape = function({ x, y, z } = halfExtents) {
  const btHalfExtents = new Ammo.btVector3(x, y, z);
  const collisionShape = new Ammo.btBoxShape(btHalfExtents);
  Ammo.destroy(btHalfExtents);
  return collisionShape;
};

const _createCylinderShape = function({ x, y, z } = halfExtents, cylinderAxis) {
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

const _createConeShape = function({ x, y, z } = halfExtents, cylinderAxis) {
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

const _createCapsuleShape = function({ x, y, z } = halfExtents, capsuleAxis) {
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
  const vertex = new THREE.Vector3();
  return function(root, mergeGeometry, margin, maxVertices) {
    const btVertex = new Ammo.btVector3();
    const originalHull = new Ammo.btConvexHullShape();
    originalHull.setMargin(margin);

    let vertexCount = 0;
    _iterateGeometries(root, mergeGeometry, geo => {
      vertexCount += geo.attributes.position.array.length / 3;
    });

    // todo: might want to implement this in a deterministic way that doesn't do O(verts) calls to Math.random
    if (vertexCount > maxVertices) {
      console.warn(`too many vertices for hull shape; sampling ~${maxVertices} from ~${vertexCount} vertices`);
    }
    const p = Math.min(1, maxVertices / vertexCount);

    _iterateGeometries(root, mergeGeometry, (geo, transform) => {
      const components = geo.attributes.position.array;
      for (let i = 0; i < components.length; i += 3) {
        if (Math.random() <= p) {
          vertex.set(components[i], components[i + 1], components[i + 2]).applyMatrix4(transform);
          btVertex.setValue(vertex.x, vertex.y, vertex.z);
          originalHull.addPoint(btVertex, i === components.length - 3); // todo: better to recalc AABB only on last geometry
        }
      }
    });

    let collisionShape = originalHull;
    if (originalHull.getNumVertices() >= 100) {
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

    root.matrixWorld.decompose(pos, quat, scale);
    const localScale = new Ammo.btVector3(scale.x, scale.y, scale.z);
    collisionShape.setLocalScaling(localScale);

    Ammo.destroy(localScale);
    Ammo.destroy(btVertex);
    return collisionShape;
  };
})();

const _createTriMeshShape = (function() {
  const pos = new THREE.Vector3();
  const quat = new THREE.Quaternion();
  const scale = new THREE.Vector3();
  const va = new THREE.Vector3();
  const vb = new THREE.Vector3();
  const vc = new THREE.Vector3();
  return function(root, mergeGeometry) {
    // todo: limit number of triangles?
    root.matrixWorld.decompose(pos, quat, scale);
    const localScale = new Ammo.btVector3(scale.x, scale.y, scale.z);
    const bta = new Ammo.btVector3();
    const btb = new Ammo.btVector3();
    const btc = new Ammo.btVector3();
    const triMesh = new Ammo.btTriangleMesh(true, true);
    triMesh.setScaling(localScale);

    _iterateGeometries(root, mergeGeometry, (geo, transform) => {
      const components = geo.attributes.position.array;
      if (geo.index) {
        for (let i = 0; i < geo.index.length; i += 3) {
          const ai = geo.index[i];
          const bi = geo.index[i + 1];
          const ci = geo.index[i + 2];
          va.set(components[ai], components[ai + 1], components[ai + 2]).applyMatrix4(transform);
          vb.set(components[bi], components[bi + 1], components[bi + 2]).applyMatrix4(transform);
          vc.set(components[ci], components[ci + 1], components[ci + 2]).applyMatrix4(transform);
          bta.setValue(va.x, va.y, va.z);
          btb.setValue(vb.x, vb.y, vb.z);
          btc.setValue(vc.x, vc.y, vc.z);
          triMesh.addTriangle(bta, btb, btc, false);
        }
      } else {
        for (let i = 0; i < components.length; i += 9) {
          va.set(components[i + 0], components[i + 1], components[i + 2]).applyMatrix4(transform);
          vb.set(components[i + 3], components[i + 4], components[i + 5]).applyMatrix4(transform);
          vc.set(components[i + 6], components[i + 7], components[i + 8]).applyMatrix4(transform);
          bta.setValue(va.x, va.y, va.z);
          btb.setValue(vb.x, vb.y, vb.z);
          btc.setValue(vc.x, vc.y, vc.z);
          triMesh.addTriangle(bta, btb, btc, false);
        }
      }
    });

    const collisionShape = new Ammo.btBvhTriangleMeshShape(triMesh, true, true);
    collisionShape.resources = [triMesh];

    Ammo.destroy(localScale);
    Ammo.destroy(bta);
    Ammo.destroy(btb);
    Ammo.destroy(btc);
    return collisionShape;
  };
})();
