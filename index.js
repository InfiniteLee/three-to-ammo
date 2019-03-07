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

    let collisionShape;
    let meshes;

    if ((mergeGeometry || autoGenerateShape) && !sceneRoot) {
      console.warn("cannot use mergeGeometry, autoGenerateShape, or recenter if sceneRoot is null");
      return null;
    }

    if (mergeGeometry) {
      meshes = _getMeshes(sceneRoot);
    } else {
      meshes = [sceneRoot];
    }

    const computeRadius = function() {
      computeBounds(sceneRoot, meshes, bounds);
      computeSphere(sceneRoot, meshes, bounds, sphere);
      return sphere.radius;
    };

    const computeHalfExtents = function() {
      computeBounds(sceneRoot, meshes, bounds);
      return halfExtents.subVectors(bounds.max, bounds.min).multiplyScalar(0.5).clampScalar(minHalfExtent, maxHalfExtent);
    };

    //TODO: Support convex hull decomposition, compound shapes, gimpact (dynamic trimesh)
    switch (type) {
    case Type.BOX: {
      const hx = options.halfExtents || computeHalfExtents();
      collisionShape = _createBoxShape(hx);
      break;
    }
    case Type.CYLINDER: {
      const hx = options.halfExtents || computeHalfExtents();
      collisionShape = _createCylinderShape(hx, cylinderAxis);
      break;
    }
    case Type.CAPSULE: {
      const hx = options.halfExtents || computeHalfExtents();
      collisionShape = _createCapsuleShape(hx, cylinderAxis);
      break;
    }
    case Type.CONE: {
      const hx = options.halfExtents || computeHalfExtents();
      collisionShape = _createConeShape(hx, cylinderAxis);
      break;
    }
    case Type.SPHERE: {
      const radius = !isNaN(options.sphereRadius) ? options.sphereRadius : computeRadius();
      collisionShape = new Ammo.btSphereShape(radius);
      break;
    }
    case Type.HULL: {
      collisionShape = _createHullShape(sceneRoot, meshes, margin, options.hullMaxVertices || 100000);
      break;
    }
    case Type.MESH: {
      collisionShape = _createTriMeshShape(sceneRoot, meshes);
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
  const matrix = new THREE.Matrix4();
  const inverse = new THREE.Matrix4();
  const vertex = new THREE.Vector3();
  return function(sceneRoot, meshes, margin, maxVertices) {
    sceneRoot.matrixWorld.decompose(pos, quat, scale);
    const localScale = new Ammo.btVector3(scale.x, scale.y, scale.z);
    const btVertex = new Ammo.btVector3();
    const originalHull = new Ammo.btConvexHullShape();
    originalHull.setMargin(margin);
    inverse.getInverse(sceneRoot.matrixWorld);

    let vertexCount = 0;
    for (let mesh of meshes) {
      const geometry = getBufferGeometry(mesh.geometry);
      const components = geometry.attributes.position.array;
      vertexCount += components.length / 3;
    }

    console.log(`Creating hull shape with ${vertexCount} vertices.`);
    if (vertexCount > maxVertices) {
      console.warn(`too many vertices for hull shape; sampling ~${maxVertices} from ~${vertexCount} vertices`);
    }
    const p = Math.min(1, maxVertices / vertexCount);

    for (let mi = 0; mi < meshes.length; mi++) {
      const mesh = meshes[mi];
      const isLastMesh = mi === meshes.length - 1;
      const geometry = getBufferGeometry(mesh.geometry);
      const components = geometry.attributes.position.array;
      if (mesh !== sceneRoot) {
        matrix.multiplyMatrices(inverse, mesh.matrixWorld);
      } else {
        matrix.identity();
      }
      for (let i = 0; i < components.length; i += 3) {
        if (Math.random() <= p) {
          vertex.set(components[i], components[i + 1], components[i + 2]).applyMatrix4(matrix);
          btVertex.setValue(vertex.x, vertex.y, vertex.z);
          originalHull.addPoint(btVertex, i === components.length - 3);
        }
      }
    }

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
  const matrix = new THREE.Matrix4();
  const inverse = new THREE.Matrix4();
  const va = new THREE.Vector3();
  const vb = new THREE.Vector3();
  const vc = new THREE.Vector3();
  return function(sceneRoot, meshes) {
    //TODO: limit number of triangles?
    const bta = new Ammo.btVector3();
    const btb = new Ammo.btVector3();
    const btc = new Ammo.btVector3();
    const triMesh = new Ammo.btTriangleMesh(true, false);
    inverse.getInverse(sceneRoot.matrixWorld);

    for (let mi = 0; mi < meshes.length; mi++) {
      const mesh = meshes[mi];
      const isLastMesh = mi === meshes.length - 1;
      const geometry = getBufferGeometry(mesh.geometry);
      const components = geometry.attributes.position.array;
      if (mesh !== sceneRoot) {
        matrix.multiplyMatrices(inverse, mesh.matrixWorld);
      } else {
        matrix.identity();
      }
      if (geometry.index) {
        for (let i = 0; i < geometry.index.length; i += 3) {
          const ai = geometry.index[i];
          const bi = geometry.index[i + 1];
          const ci = geometry.index[i + 2];
          va.set(components[ai], components[ai + 1], components[ai + 2]).applyMatrix4(matrix);
          vb.set(components[bi], components[bi + 1], components[bi + 2]).applyMatrix4(matrix);
          vc.set(components[ci], components[ci + 1], components[ci + 2]).applyMatrix4(matrix);
          bta.setValue(va.x, va.y, va.z);
          btb.setValue(vb.x, vb.y, vb.z);
          btc.setValue(vc.x, vc.y, vc.z);
          triMesh.addTriangle(bta, btb, btc, isLastMesh && i == geometry.index.length - 3);
        }
      } else {
        for (let i = 0; i < components.length; i += 9) {
          va.set(components[i + 0], components[i + 1], components[i + 2]).applyMatrix4(matrix);
          vb.set(components[i + 3], components[i + 4], components[i + 5]).applyMatrix4(matrix);
          vc.set(components[i + 6], components[i + 7], components[i + 8]).applyMatrix4(matrix);
          bta.setValue(va.x, va.y, va.z);
          btb.setValue(vb.x, vb.y, vb.z);
          btc.setValue(vc.x, vc.y, vc.z);
          triMesh.addTriangle(bta, btb, btc, isLastMesh && i == components.length - 9);
        }
      }
    }

    sceneRoot.matrixWorld.decompose(pos, quat, scale);
    const localScale = new Ammo.btVector3(scale.x, scale.y, scale.z);
    const collisionShape = new Ammo.btBvhTriangleMeshShape(triMesh, true, true);
    collisionShape.setLocalScaling(localScale);
    collisionShape.resources = [triMesh];
    Ammo.destroy(localScale);
    Ammo.destroy(bta);
    Ammo.destroy(btb);
    Ammo.destroy(btc);
    return collisionShape;
  };
})();

const getBufferGeometry = (function() {
  const bufferGeometry = new THREE.BufferGeometry();
  return function(geo) {
    return geo.isBufferGeometry ? geo : bufferGeometry.fromGeometry(geo);
  };
})();

const computeSphere = (function() {
  const matrix = new THREE.Matrix4();
  const inverse = new THREE.Matrix4();
  const vertex = new THREE.Vector3();
  return function(sceneRoot, meshes, bounds, target) {

    let { x: cx, y: cy, z: cz } = bounds.getCenter(target.center);
    let maxRadiusSq = 0;
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
        vertex.set(components[i], components[i + 1], components[i + 2]).applyMatrix4(matrix);
        const dx = cx - vertex.x;
        const dy = cy - vertex.y;
        const dz = cz - vertex.z;
        maxRadiusSq = Math.max(maxRadiusSq, dx * dx + dy * dy + dz * dz);
      }
    }
    target.radius = Math.sqrt(maxRadiusSq);
    return target;
  };
})();

const computeBounds = (function() {
  const v = new THREE.Vector3();
  const matrix = new THREE.Matrix4();
  const inverse = new THREE.Matrix4();
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
        v.set(components[i], components[i + 1], components[i + 2]).applyMatrix4(matrix);
		if ( v.x < minX ) minX = v.x;
		if ( v.y < minY ) minY = v.y;
		if ( v.z < minZ ) minZ = v.z;
		if ( v.x > maxX ) maxX = v.x;
		if ( v.y > maxY ) maxY = v.y;
		if ( v.z > maxZ ) maxZ = v.z;
      }
    }

    target.min.set(minX, minY, minZ);
    target.max.set(maxX, maxY, maxZ);
    return target;
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
