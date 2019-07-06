"use strict";
/* global Ammo,THREE */

const TYPE = (exports.TYPE = {
  BOX: "box",
  CYLINDER: "cylinder",
  SPHERE: "sphere",
  CAPSULE: "capsule",
  CONE: "cone",
  HULL: "hull",
  HACD: "hacd", //Hierarchical Approximate Convex Decomposition
  VHACD: "vhacd", //Volumetric Hierarchical Approximate Convex Decomposition
  MESH: "mesh",
  HEIGHTFIELD: "heightfield"
});

const FIT = (exports.FIT = {
  ALL: "all", //A single shape is automatically sized to bound all meshes within the entity.
  MANUAL: "manual" //A single shape is sized manually. Requires halfExtents or sphereRadius.
});

const HEIGHTFIELD_DATA_TYPE = (exports.HEIGHTFIELD_DATA_TYPE = {
  short: "short",
  float: "float"
});

const hasUpdateMatricesFunction = THREE.Object3D.prototype.hasOwnProperty("updateMatrices");

exports.createCollisionShapes = function(root, options) {
  switch (options.type) {
    case TYPE.BOX:
      return [this.createBoxShape(root, options)];
    case TYPE.CYLINDER:
      return [this.createCylinderShape(root, options)];
    case TYPE.CAPSULE:
      return [this.createCapsuleShape(root, options)];
    case TYPE.CONE:
      return [this.createConeShape(root, options)];
    case TYPE.SPHERE:
      return [this.createSphereShape(root, options)];
    case TYPE.HULL:
      return [this.createHullShape(root, options)];
    case TYPE.HACD:
      return this.createHACDShapes(root, options);
    case TYPE.VHACD:
      return this.createVHACDShapes(root, options);
    case TYPE.MESH:
      return [this.createTriMeshShape(root, options)];
    case TYPE.HEIGHTFIELD:
      return [this.createHeightfieldTerrainShape(root, options)];
    default:
      console.warn(options.type + " is not currently supported");
      return [];
  }
};

//TODO: support gimpact (dynamic trimesh) and heightmap

exports.createBoxShape = function(root, options) {
  options.type = TYPE.BOX;
  _setOptions(options);

  if (options.fit === FIT.ALL) {
    options.halfExtents = _computeHalfExtents(
      root,
      _computeBounds(root, options),
      options.minHalfExtent,
      options.maxHalfExtent
    );
  }

  const btHalfExtents = new Ammo.btVector3(options.halfExtents.x, options.halfExtents.y, options.halfExtents.z);
  const collisionShape = new Ammo.btBoxShape(btHalfExtents);
  Ammo.destroy(btHalfExtents);

  _finishCollisionShape(collisionShape, options, _computeScale(root, options));
  return collisionShape;
};

exports.createCylinderShape = function(root, options) {
  options.type = TYPE.CYLINDER;
  _setOptions(options);

  if (options.fit === FIT.ALL) {
    options.halfExtents = _computeHalfExtents(
      root,
      _computeBounds(root, options),
      options.minHalfExtent,
      options.maxHalfExtent
    );
  }

  const btHalfExtents = new Ammo.btVector3(options.halfExtents.x, options.halfExtents.y, options.halfExtents.z);
  const collisionShape = (() => {
    switch (options.cylinderAxis) {
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

  _finishCollisionShape(collisionShape, options, _computeScale(root, options));
  return collisionShape;
};

exports.createCapsuleShape = function(root, options) {
  options.type = TYPE.CAPSULE;
  _setOptions(options);

  if (options.fit === FIT.ALL) {
    options.halfExtents = _computeHalfExtents(
      root,
      _computeBounds(root, options),
      options.minHalfExtent,
      options.maxHalfExtent
    );
  }

  const { x, y, z } = options.halfExtents;
  const collisionShape = (() => {
    switch (options.cylinderAxis) {
      case "y":
        return new Ammo.btCapsuleShape(Math.max(x, z), y * 2);
      case "x":
        return new Ammo.btCapsuleShapeX(Math.max(y, z), x * 2);
      case "z":
        return new Ammo.btCapsuleShapeZ(Math.max(x, y), z * 2);
    }
    return null;
  })();

  _finishCollisionShape(collisionShape, options, _computeScale(root, options));
  return collisionShape;
};

exports.createConeShape = function(root, options) {
  options.type = TYPE.CONE;
  _setOptions(options);

  if (options.fit === FIT.ALL) {
    options.halfExtents = _computeHalfExtents(
      root,
      _computeBounds(root, options),
      options.minHalfExtent,
      options.maxHalfExtent
    );
  }

  const { x, y, z } = options.halfExtents;
  const collisionShape = (() => {
    switch (options.cylinderAxis) {
      case "y":
        return new Ammo.btConeShape(Math.max(x, z), y * 2);
      case "x":
        return new Ammo.btConeShapeX(Math.max(y, z), x * 2);
      case "z":
        return new Ammo.btConeShapeZ(Math.max(x, y), z * 2);
    }
    return null;
  })();

  _finishCollisionShape(collisionShape, options, _computeScale(root, options));
  return collisionShape;
};

exports.createSphereShape = function(root, options) {
  options.type = TYPE.SPHERE;
  _setOptions(options);

  let radius;
  if (options.fit === FIT.MANUAL && !isNaN(options.sphereRadius)) {
    radius = options.sphereRadius;
  } else {
    radius = _computeRadius(root, options, _computeBounds(root, options));
  }

  const collisionShape = new Ammo.btSphereShape(radius);
  _finishCollisionShape(collisionShape, options, _computeScale(root, options));

  return collisionShape;
};

exports.createHullShape = (function() {
  const vertex = new THREE.Vector3();
  const center = new THREE.Vector3();
  return function(root, options) {
    options.type = TYPE.HULL;
    _setOptions(options);

    if (options.fit === FIT.MANUAL) {
      console.warn("cannot use fit: manual with type: hull");
      return null;
    }

    const bounds = _computeBounds(root, options);

    const btVertex = new Ammo.btVector3();
    const originalHull = new Ammo.btConvexHullShape();
    originalHull.setMargin(options.margin);
    center.addVectors(bounds.max, bounds.min).multiplyScalar(0.5);

    let vertexCount = 0;
    _iterateGeometries(root, options, geo => {
      vertexCount += geo.attributes.position.array.length / 3;
    });

    const maxVertices = options.hullMaxVertices || 100000;
    // todo: might want to implement this in a deterministic way that doesn't do O(verts) calls to Math.random
    if (vertexCount > maxVertices) {
      console.warn(`too many vertices for hull shape; sampling ~${maxVertices} from ~${vertexCount} vertices`);
    }
    const p = Math.min(1, maxVertices / vertexCount);

    _iterateGeometries(root, options, (geo, transform) => {
      const components = geo.attributes.position.array;
      for (let i = 0; i < components.length; i += 3) {
        if (Math.random() <= p) {
          vertex
            .set(components[i], components[i + 1], components[i + 2])
            .applyMatrix4(transform)
            .sub(center);
          btVertex.setValue(vertex.x, vertex.y, vertex.z);
          originalHull.addPoint(btVertex, i === components.length - 3); // todo: better to recalc AABB only on last geometry
        }
      }
    });

    let collisionShape = originalHull;
    if (originalHull.getNumVertices() >= 100) {
      //Bullet documentation says don't use convexHulls with 100 verts or more
      const shapeHull = new Ammo.btShapeHull(originalHull);
      shapeHull.buildHull(options.margin);
      Ammo.destroy(originalHull);
      collisionShape = new Ammo.btConvexHullShape(
        Ammo.getPointer(shapeHull.getVertexPointer()),
        shapeHull.numVertices()
      );
      Ammo.destroy(shapeHull); // btConvexHullShape makes a copy
    }

    Ammo.destroy(btVertex);

    _finishCollisionShape(collisionShape, options, _computeScale(root, options));
    return collisionShape;
  };
})();

exports.createHACDShapes = (function() {
  const v = new THREE.Vector3();
  const center = new THREE.Vector3();
  return function(root, options) {
    options.type = TYPE.HACD;
    _setOptions(options);

    if (options.fit === FIT.MANUAL) {
      console.warn("cannot use fit: manual with type: hacd");
      return [];
    }

    if (!Ammo.hasOwnProperty("HACD")) {
      console.warn(
        "HACD unavailable in included build of Ammo.js. Visit https://github.com/mozillareality/ammo.js for the latest version."
      );
      return [];
    }

    const bounds = _computeBounds(root);
    const scale = _computeScale(root, options);

    let vertexCount = 0;
    let triCount = 0;
    center.addVectors(bounds.max, bounds.min).multiplyScalar(0.5);

    _iterateGeometries(root, options, geo => {
      vertexCount += geo.attributes.position.array.length / 3;
      if (geo.index) {
        triCount += geo.index.array.length / 3;
      } else {
        triCount += geo.attributes.position.array.length / 9;
      }
    });

    const hacd = new Ammo.HACD();
    if (options.hasOwnProperty("compacityWeight")) hacd.SetCompacityWeight(options.compacityWeight);
    if (options.hasOwnProperty("volumeWeight")) hacd.SetVolumeWeight(options.volumeWeight);
    if (options.hasOwnProperty("nClusters")) hacd.SetNClusters(options.nClusters);
    if (options.hasOwnProperty("nVerticesPerCH")) hacd.SetNVerticesPerCH(options.nVerticesPerCH);
    if (options.hasOwnProperty("concavity")) hacd.SetConcavity(options.concavity);

    const points = Ammo._malloc(vertexCount * 3 * 8);
    const triangles = Ammo._malloc(triCount * 3 * 4);
    hacd.SetPoints(points);
    hacd.SetTriangles(triangles);
    hacd.SetNPoints(vertexCount);
    hacd.SetNTriangles(triCount);

    const pptr = points / 8,
      tptr = triangles / 4;
    _iterateGeometries(root, options, (geo, transform) => {
      const components = geo.attributes.position.array;
      const indices = geo.index ? geo.index.array : null;
      for (let i = 0; i < components.length; i += 3) {
        v.set(components[i + 0], components[i + 1], components[i + 2])
          .applyMatrix4(transform)
          .sub(center);
        Ammo.HEAPF64[pptr + i + 0] = v.x;
        Ammo.HEAPF64[pptr + i + 1] = v.y;
        Ammo.HEAPF64[pptr + i + 2] = v.z;
      }
      if (indices) {
        for (let i = 0; i < indices.length; i++) {
          Ammo.HEAP32[tptr + i] = indices[i];
        }
      } else {
        for (let i = 0; i < components.length / 3; i++) {
          Ammo.HEAP32[tptr + i] = i;
        }
      }
    });

    hacd.Compute();
    Ammo._free(points);
    Ammo._free(triangles);
    const nClusters = hacd.GetNClusters();

    const shapes = [];
    for (let i = 0; i < nClusters; i++) {
      const hull = new Ammo.btConvexHullShape();
      hull.setMargin(options.margin);
      const nPoints = hacd.GetNPointsCH(i);
      const nTriangles = hacd.GetNTrianglesCH(i);
      const hullPoints = Ammo._malloc(nPoints * 3 * 8);
      const hullTriangles = Ammo._malloc(nTriangles * 3 * 4);
      hacd.GetCH(i, hullPoints, hullTriangles);

      const pptr = hullPoints / 8;
      for (let pi = 0; pi < nPoints; pi++) {
        const btVertex = new Ammo.btVector3();
        const px = Ammo.HEAPF64[pptr + pi * 3 + 0];
        const py = Ammo.HEAPF64[pptr + pi * 3 + 1];
        const pz = Ammo.HEAPF64[pptr + pi * 3 + 2];
        btVertex.setValue(px, py, pz);
        hull.addPoint(btVertex, pi === nPoints - 1);
        Ammo.destroy(btVertex);
      }

      _finishCollisionShape(hull, options, scale);
      shapes.push(hull);
    }

    return shapes;
  };
})();

exports.createVHACDShapes = (function() {
  const v = new THREE.Vector3();
  const center = new THREE.Vector3();
  return function(root, options) {
    options.type = TYPE.VHACD;
    _setOptions(options);

    if (options.fit === FIT.MANUAL) {
      console.warn("cannot use fit: manual with type: vhacd");
      return [];
    }

    if (!Ammo.hasOwnProperty("VHACD")) {
      console.warn(
        "VHACD unavailable in included build of Ammo.js. Visit https://github.com/mozillareality/ammo.js for the latest version."
      );
      return [];
    }

    const bounds = _computeBounds(root, options);
    const scale = _computeScale(root, options);

    let vertexCount = 0;
    let triCount = 0;
    center.addVectors(bounds.max, bounds.min).multiplyScalar(0.5);

    _iterateGeometries(root, options, geo => {
      vertexCount += geo.attributes.position.count;
      if (geo.index) {
        triCount += geo.index.count / 3;
      } else {
        triCount += geo.attributes.position.count / 3;
      }
    });

    const vhacd = new Ammo.VHACD();
    const params = new Ammo.Parameters();
    //https://kmamou.blogspot.com/2014/12/v-hacd-20-parameters-description.html
    if (options.hasOwnProperty("resolution")) params.set_m_resolution(options.resolution);
    if (options.hasOwnProperty("depth")) params.set_m_depth(options.depth);
    if (options.hasOwnProperty("concavity")) params.set_m_concavity(options.concavity);
    if (options.hasOwnProperty("planeDownsampling")) params.set_m_planeDownsampling(options.planeDownsampling);
    if (options.hasOwnProperty("convexhullDownsampling"))
      params.set_m_convexhullDownsampling(options.convexhullDownsampling);
    if (options.hasOwnProperty("alpha")) params.set_m_alpha(options.alpha);
    if (options.hasOwnProperty("beta")) params.set_m_beta(options.beta);
    if (options.hasOwnProperty("gamma")) params.set_m_gamma(options.gamma);
    if (options.hasOwnProperty("pca")) params.set_m_pca(options.pca);
    if (options.hasOwnProperty("mode")) params.set_m_mode(options.mode);
    if (options.hasOwnProperty("maxNumVerticesPerCH")) params.set_m_maxNumVerticesPerCH(options.maxNumVerticesPerCH);
    if (options.hasOwnProperty("minVolumePerCH")) params.set_m_minVolumePerCH(options.minVolumePerCH);
    if (options.hasOwnProperty("convexhullApproximation"))
      params.set_m_convexhullApproximation(options.convexhullApproximation);
    if (options.hasOwnProperty("oclAcceleration")) params.set_m_oclAcceleration(options.oclAcceleration);

    const points = Ammo._malloc(vertexCount * 3 * 8);
    const triangles = Ammo._malloc(triCount * 3 * 4);

    let pptr = points / 8,
      tptr = triangles / 4;
    _iterateGeometries(root, options, (geo, transform) => {
      const components = geo.attributes.position.array;
      const indices = geo.index ? geo.index.array : null;
      for (let i = 0; i < components.length; i += 3) {
        v.set(components[i + 0], components[i + 1], components[i + 2])
          .applyMatrix4(transform)
          .sub(center);
        Ammo.HEAPF64[pptr + 0] = v.x;
        Ammo.HEAPF64[pptr + 1] = v.y;
        Ammo.HEAPF64[pptr + 2] = v.z;
        pptr += 3;
      }
      if (indices) {
        for (let i = 0; i < indices.length; i++) {
          Ammo.HEAP32[tptr] = indices[i];
          tptr++;
        }
      } else {
        for (let i = 0; i < components.length / 3; i++) {
          Ammo.HEAP32[tptr] = i;
          tptr++;
        }
      }
    });

    vhacd.Compute(points, 3, vertexCount, triangles, 3, triCount, params);
    Ammo._free(points);
    Ammo._free(triangles);
    const nHulls = vhacd.GetNConvexHulls();

    const shapes = [];
    const ch = new Ammo.ConvexHull();
    for (let i = 0; i < nHulls; i++) {
      vhacd.GetConvexHull(i, ch);
      const nPoints = ch.get_m_nPoints();
      const hullPoints = ch.get_m_points();

      const hull = new Ammo.btConvexHullShape();
      hull.setMargin(options.margin);

      for (let pi = 0; pi < nPoints; pi++) {
        const btVertex = new Ammo.btVector3();
        const px = ch.get_m_points(pi * 3 + 0);
        const py = ch.get_m_points(pi * 3 + 1);
        const pz = ch.get_m_points(pi * 3 + 2);
        btVertex.setValue(px, py, pz);
        hull.addPoint(btVertex, pi === nPoints - 1);
        Ammo.destroy(btVertex);
      }

      _finishCollisionShape(hull, options, scale);
      shapes.push(hull);
    }
    Ammo.destroy(ch);
    Ammo.destroy(vhacd);

    return shapes;
  };
})();

exports.createTriMeshShape = (function() {
  const va = new THREE.Vector3();
  const vb = new THREE.Vector3();
  const vc = new THREE.Vector3();
  return function(root, options) {
    options.type = TYPE.MESH;
    _setOptions(options);

    if (options.fit === FIT.MANUAL) {
      console.warn("cannot use fit: manual with type: mesh");
      return null;
    }

    const scale = _computeScale(root, options);

    const bta = new Ammo.btVector3();
    const btb = new Ammo.btVector3();
    const btc = new Ammo.btVector3();
    const triMesh = new Ammo.btTriangleMesh(true, false);

    _iterateGeometries(root, options, (geo, transform) => {
      const components = geo.attributes.position.array;
      if (geo.index) {
        for (let i = 0; i < geo.index.count; i += 3) {
          const ai = geo.index.array[i] * 3;
          const bi = geo.index.array[i + 1] * 3;
          const ci = geo.index.array[i + 2] * 3;
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

    const localScale = new Ammo.btVector3(scale.x, scale.y, scale.z);
    triMesh.setScaling(localScale);
    Ammo.destroy(localScale);

    const collisionShape = new Ammo.btBvhTriangleMeshShape(triMesh, true, true);
    collisionShape.resources = [triMesh];

    Ammo.destroy(bta);
    Ammo.destroy(btb);
    Ammo.destroy(btc);

    _finishCollisionShape(collisionShape, options);
    return collisionShape;
  };
})();

exports.createHeightfieldTerrainShape = function(root, options) {
  _setOptions(options);

  if (options.fit === FIT.ALL) {
    console.warn("cannot use fit: all with type: heightfield");
    return null;
  }
  const heightfieldDistance = options.heightfieldDistance || 1;
  const heightfieldData = options.heightfieldData || [];
  const heightScale = options.heightScale || 0;
  const upAxis = options.hasOwnProperty("upAxis") ? options.upAxis : 1; // x = 0; y = 1; z = 2
  const hdt = (() => {
    switch (options.heightDataType) {
      case "short":
        return Ammo.PHY_SHORT;
      case "float":
        return Ammo.PHY_FLOAT;
      default:
        return Ammo.PHY_FLOAT;
    }
  })();
  const flipQuadEdges = options.hasOwnProperty("flipQuadEdges") ? options.flipQuadEdges : true;

  const heightStickLength = heightfieldData.length;
  const heightStickWidth = heightStickLength > 0 ? heightfieldData[0].length : 0;

  const data = Ammo._malloc(heightStickLength * heightStickWidth * 4);
  const ptr = data / 4;

  let minHeight = Number.POSITIVE_INFINITY;
  let maxHeight = Number.NEGATIVE_INFINITY;
  let index = 0;
  for (let l = 0; l < heightStickLength; l++) {
    for (let w = 0; w < heightStickWidth; w++) {
      const height = heightfieldData[l][w];
      Ammo.HEAPF32[ptr + index] = height;
      index++;
      minHeight = Math.min(minHeight, height);
      maxHeight = Math.max(maxHeight, height);
    }
  }

  const collisionShape = new Ammo.btHeightfieldTerrainShape(
    heightStickWidth,
    heightStickLength,
    data,
    heightScale,
    minHeight,
    maxHeight,
    upAxis,
    hdt,
    flipQuadEdges
  );

  const scale = new Ammo.btVector3(heightfieldDistance, 1, heightfieldDistance);
  collisionShape.setLocalScaling(scale);
  Ammo.destroy(scale);

  collisionShape.heightfieldData = data;

  _finishCollisionShape(collisionShape, options);
  return collisionShape;
};

function _setOptions(options) {
  options.fit = options.hasOwnProperty("fit") ? options.fit : FIT.ALL;
  options.type = options.type || TYPE.HULL;
  options.minHalfExtent = options.hasOwnProperty("minHalfExtent") ? options.minHalfExtent : 0;
  options.maxHalfExtent = options.hasOwnProperty("maxHalfExtent") ? options.maxHalfExtent : Number.POSITIVE_INFINITY;
  options.cylinderAxis = options.cylinderAxis || "y";
  options.margin = options.hasOwnProperty("margin") ? options.margin : 0.01;
  options.includeInvisible = options.hasOwnProperty("includeInvisible") ? options.includeInvisible : false;

  if (!options.offset) {
    options.offset = new THREE.Vector3();
  }

  if (!options.orientation) {
    options.orientation = new THREE.Quaternion();
  }
}

const _finishCollisionShape = function(collisionShape, options, scale) {
  collisionShape.type = options.type;
  collisionShape.setMargin(options.margin);
  collisionShape.destroy = () => {
    for (let res of collisionShape.resources || []) {
      Ammo.destroy(res);
    }
    if (collisionShape.heightfieldData) {
      Ammo._free(collisionShape.heightfieldData);
    }
    Ammo.destroy(collisionShape);
  };

  const localTransform = new Ammo.btTransform();
  const rotation = new Ammo.btQuaternion();
  localTransform.setIdentity();

  localTransform.getOrigin().setValue(options.offset.x, options.offset.y, options.offset.z);
  rotation.setValue(options.orientation.x, options.orientation.y, options.orientation.z, options.orientation.w);

  localTransform.setRotation(rotation);
  Ammo.destroy(rotation);

  if (scale) {
    const localScale = new Ammo.btVector3(scale.x, scale.y, scale.z);
    collisionShape.setLocalScaling(localScale);
    Ammo.destroy(localScale);
  }

  collisionShape.localTransform = localTransform;
};

// Calls `cb(geo, transform)` for each geometry under `root` whose vertices we should take into account for the physics shape.
// `transform` is the transform required to transform the given geometry's vertices into root-local space.
const _iterateGeometries = (function() {
  const transform = new THREE.Matrix4();
  const inverse = new THREE.Matrix4();
  const bufferGeometry = new THREE.BufferGeometry();
  return function(root, options, cb) {
    inverse.getInverse(root.matrixWorld);
    root.traverse(mesh => {
      if (
        mesh.isMesh &&
        (!THREE.Sky || mesh.__proto__ != THREE.Sky.prototype) &&
        (options.includeInvisible || ((mesh.el && mesh.el.object3D.visible) || mesh.visible))
      ) {
        if (mesh === root) {
          transform.identity();
        } else {
          if (hasUpdateMatricesFunction) mesh.updateMatrices();
          transform.multiplyMatrices(inverse, mesh.matrixWorld);
        }
        // todo: might want to return null xform if this is the root so that callers can avoid multiplying
        // things by the identity matrix
        cb(mesh.geometry.isBufferGeometry ? mesh.geometry : bufferGeometry.fromGeometry(mesh.geometry), transform);
      }
    });
  };
})();

const _computeScale = function(root, options) {
  const scale = new THREE.Vector3(1, 1, 1);
  if (options.fit === FIT.ALL) {
    scale.setFromMatrixScale(root.matrixWorld);
  }
  return scale;
};

const _computeRadius = (function() {
  const v = new THREE.Vector3();
  const center = new THREE.Vector3();
  return function(root, options, bounds) {
    let maxRadiusSq = 0;
    let { x: cx, y: cy, z: cz } = bounds.getCenter(center);
    _iterateGeometries(root, options, (geo, transform) => {
      const components = geo.attributes.position.array;
      for (let i = 0; i < components.length; i += 3) {
        v.set(components[i], components[i + 1], components[i + 2]).applyMatrix4(transform);
        const dx = cx - v.x;
        const dy = cy - v.y;
        const dz = cz - v.z;
        maxRadiusSq = Math.max(maxRadiusSq, dx * dx + dy * dy + dz * dz);
      }
    });
    return Math.sqrt(maxRadiusSq);
  };
})();

const _computeHalfExtents = function(root, bounds, minHalfExtent, maxHalfExtent) {
  const halfExtents = new THREE.Vector3();
  return halfExtents
    .subVectors(bounds.max, bounds.min)
    .multiplyScalar(0.5)
    .clampScalar(minHalfExtent, maxHalfExtent);
};

const _computeLocalOffset = function(matrix, bounds, target) {
  target
    .addVectors(bounds.max, bounds.min)
    .multiplyScalar(0.5)
    .applyMatrix4(matrix);
  return target;
};

// returns the bounding box for the geometries underneath `root`.
const _computeBounds = (function() {
  const v = new THREE.Vector3();
  return function(root, options) {
    const bounds = new THREE.Box3();
    let minX = +Infinity;
    let minY = +Infinity;
    let minZ = +Infinity;
    let maxX = -Infinity;
    let maxY = -Infinity;
    let maxZ = -Infinity;
    bounds.min.set(0, 0, 0);
    bounds.max.set(0, 0, 0);
    _iterateGeometries(root, options, (geo, transform) => {
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
    bounds.min.set(minX, minY, minZ);
    bounds.max.set(maxX, maxY, maxZ);
    return bounds;
  };
})();
