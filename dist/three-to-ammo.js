"use strict";

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

/* global Ammo,THREE */

var SHAPE_BOX = "box";
var SHAPE_CYLINDER = "cylinder";
var SHAPE_SPHERE = "sphere";
var SHAPE_CAPSULE = "capsule";
var SHAPE_CONE = "cone";
var SHAPE_HULL = "hull";
var SHAPE_MESH = "mesh";

var ThreeToAmmo = function () {
  function ThreeToAmmo() {
    _classCallCheck(this, ThreeToAmmo);

    this.vertexPool = [];
    this.vertices = [];
    this.matrix = new THREE.Matrix4();
    this.inverse = new THREE.Matrix4();
    this.halfExtents = new THREE.Vector3();
    this.offset = new THREE.Matrix4();
    this.center = new THREE.Vector3();
    this.geometries = [];
    this.boundingBox = {
      min: new THREE.Vector3(Number.MAX_VALUE),
      max: new THREE.Vector3(Number.MIN_VALUE)
    };
    this.pos = new THREE.Vector3();
    this.quat = new THREE.Quaternion();
    this.boundingBox = new THREE.Box3();
  }

  _createClass(ThreeToAmmo, [{
    key: "createCollisionShape",
    value: function createCollisionShape(sceneRoot, options) {
      var autoGenerateShape = options.autoGenerateShape || true;
      var shape = options.shape || SHAPE_HULL;
      var recenter = options.recenter || false;
      if (options.halfExtents) {
        this.halfExtents.set(options.halfExtents.x, options.halfExtents.y, options.halfExtents.z);
      }
      var cylinderAxis = options.cylinderAxis || "y";
      var sphereRadius = options.sphereRadius || 1;
      var collisionShape = void 0;
      var triMesh = void 0;
      var shapeHull = void 0;

      var meshes = this._getMeshes(sceneRoot);

      if (shape !== "mesh") {
        if (recenter) {
          this._recenter(sceneRoot, meshes);
        }
      }

      var vertices = this._getVertices(sceneRoot, meshes);
      boundingBox.setFromPoints(vertices);

      //TODO: Support convex hull decomposition, compound shapes, gimpact (dynamic trimesh)
      switch (shape) {
        case "box":
          {
            var _getHalfExtents2 = this._getHalfExtents(boundingBox),
                x = _getHalfExtents2.x,
                y = _getHalfExtents2.y,
                z = _getHalfExtents2.z;

            var _halfExtents = new Ammo.btVector3(x, y, z);
            collisionShape = new Ammo.btBoxShape(_halfExtents);
            Ammo.destroy(_halfExtents);
            break;
          }
        case "sphere":
          {
            var radius = 1;
            if (sphereRadius) {
              radius = sphereRadius;
            } else {
              var sphere = new THREE.Sphere();
              sphere.setFromPoints(vertices);
              if (isFinite(sphere.radius)) {
                radius = sphere.radius;
              }
            }
            collisionShape = new Ammo.btSphereShape(radius);
            break;
          }
        case "cylinder":
          {
            var _getHalfExtents3 = this._getHalfExtents(boundingBox),
                _x = _getHalfExtents3.x,
                _y = _getHalfExtents3.y,
                _z = _getHalfExtents3.z;

            var _halfExtents2 = new Ammo.btVector3(_x, _y, _z);
            switch (cylinderAxis) {
              case "y":
                collisionShape = new Ammo.btCylinderShape(_halfExtents2);
                break;
              case "x":
                collisionShape = new Ammo.btCylinderShapeX(_halfExtents2);
                break;
              case "z":
                collisionShape = new Ammo.btCylinderShapeZ(_halfExtents2);
                break;
            }
            Ammo.destroy(_halfExtents2);
            break;
          }
        case "capsule":
          {
            var _getHalfExtents4 = this._getHalfExtents(boundingBox),
                _x2 = _getHalfExtents4.x,
                _y2 = _getHalfExtents4.y,
                _z2 = _getHalfExtents4.z;

            switch (cylinderAxis) {
              case "y":
                collisionShape = new Ammo.btCapsuleShape(Math.max(_x2, _z2), _y2 * 2);
                break;
              case "x":
                collisionShape = new Ammo.btCapsuleShapeX(Math.max(_y2, _z2), _x2 * 2);
                break;
              case "z":
                collisionShape = new Ammo.btCapsuleShapeZ(Math.max(_x2, _y2), _z2 * 2);
                break;
            }
            break;
          }
        case "cone":
          {
            var _getHalfExtents5 = this._getHalfExtents(boundingBox),
                _x3 = _getHalfExtents5.x,
                _y3 = _getHalfExtents5.y,
                _z3 = _getHalfExtents5.z;

            switch (cylinderAxis) {
              case "y":
                collisionShape = new Ammo.btConeShape(Math.max(_x3, _z3), _y3 * 2);
                break;
              case "x":
                collisionShape = new Ammo.btConeShapeX(Math.max(_y3, _z3), _x3 * 2);
                break;
              case "z":
                collisionShape = new Ammo.btConeShapeZ(Math.max(_x3, _y3), _z3 * 2);
                break;
            }
            break;
          }
        case "hull":
          {
            var scale = new Ammo.btVector3(sceneRoot.scale.x, sceneRoot.scale.y, sceneRoot.scale.z);
            var vec3 = new Ammo.btVector3();
            var originalHull = new Ammo.btConvexHullShape();
            originalHull.setMargin(data.margin);

            for (var i = 0; i < vertices.length; i++) {
              vec3.setValue(vertices[i].x, vertices[i].y, vertices[i].z);
              originalHull.addPoint(vec3, i == vertices.length - 1);
            }

            collisionShape = originalHull;
            if (originalHull.getNumVertices() >= 100) {
              //Bullet documentation says don't use convexHulls with 100 verts or more
              shapeHull = new Ammo.btShapeHull(originalHull);
              shapeHull.buildHull(data.margin);
              Ammo.destroy(originalHull);
              collisionShape = new Ammo.btConvexHullShape(Ammo.getPointer(shapeHull.getVertexPointer()), shapeHull.numVertices());
            }
            collisionShape.setLocalScaling(scale);

            Ammo.destroy(scale);
            Ammo.destroy(vec3);
            break;
          }
        case "mesh":
          {
            if (data.type !== "static") {
              //TODO: support btTriangleMeshShape for dynamic trimeshes. (not performant)
              console.warn("non-static mesh colliders are not currently supported");
              break;
            }

            var a = new Ammo.btVector3();
            var b = new Ammo.btVector3();
            var c = new Ammo.btVector3();
            triMesh = new Ammo.btTriangleMesh(true, false);

            for (var j = 0; j < vertices.length; j += 3) {
              a.setValue(vertices[j].x, vertices[j].y, vertices[j].z);
              b.setValue(vertices[j + 1].x, vertices[j + 1].y, vertices[j + 1].z);
              c.setValue(vertices[j + 2].x, vertices[j + 2].y, vertices[j + 2].z);
              triMesh.addTriangle(a, b, c, j == vertices.length - 3);
            }

            collisionShape = new Ammo.btBvhTriangleMeshShape(triMesh, true, true);
            collisionShape.setMargin(data.margin);
            //TODO: support btScaledBvhTriangleMeshShape?

            Ammo.destroy(a);
            Ammo.destroy(b);
            Ammo.destroy(c);
            break;
          }

        default:
          console.warn(data.shape + " is not currently supported");
          return;
      }

      return {
        collisionShape: collisionShape,
        destroy: function destroy() {
          if (shapeHull) Ammo.destroy(shapeHull);
          if (triMesh) Ammo.destroy(triMesh);
          if (collisionShape) Ammo.destroy(collisionShape);
        }
      };
    }
  }, {
    key: "_getMeshes",
    value: function _getMeshes(sceneRoot) {
      var meshes = [];
      sceneRoot.traverse(function (o) {
        if (o.type === "Mesh" && (!THREE.Sky || o.__proto__ != THREE.Sky.prototype)) {
          meshes.push(o);
        }
      });
      return meshes;
    }
  }, {
    key: "_getVertices",
    value: function _getVertices(sceneRoot, meshes) {
      while (this.vertices.length > 0) {
        this.vertexPool.push(vertices.pop());
      }

      this.inverse.getInverse(sceneRoot.matrixWorld);

      for (var j = 0; j < meshes.length; j++) {
        var mesh = meshes[j];

        var geometry = mesh.geometry.index ? mesh.geometry.toNonIndexed() : mesh.geometry.clone();

        if (this.shape === SHAPE_MESH) {
          geometry.applyMatrix(mesh.matrixWorld);
        } else {
          this.matrix.multiplyMatrices(this.inverse, mesh.matrixWorld);
          geometry.applyMatrix(this.matrix);
        }

        if (geometry.isBufferGeometry) {
          var components = geometry.attributes.position.array;
          for (var i = 0; i < components.length; i += 3) {
            var x = components[i];
            var y = components[i + 1];
            var z = components[i + 2];

            if (this.vertexPool.length > 0) {
              this.vertices.push(this.vertexPool.pop().set(x, y, z));
            } else {
              this.vertices.push(new THREE.Vector3(x, y, z));
            }
          }
        } else {
          for (var _i = 0; _i < geometry.vertices.length; _i++) {
            var vertex = geometry.vertices[_i];
            if (this.vertexPool.length > 0) {
              this.vertices.push(this.vertexPool.pop().copy(vertex));
            } else {
              this.vertices.push(new THREE.Vector3(vertex.x, vertex.y, vertex.z));
            }
          }
        }
      }

      return vertices;
    }
  }, {
    key: "_getHalfExtents",
    value: function _getHalfExtents(boundingBox) {
      if (this.autoGenerateShape) {
        var min = boundingBox.min,
            max = boundingBox.max;

        halfExtents.subVectors(max, min).multiplyScalar(0.5);
      }

      return {
        x: this.halfExtents.x,
        y: this.halfExtents.y,
        z: this.halfExtents.z
      };
    }
  }, {
    key: "_recenter",
    value: function _recenter(sceneRoot, meshes) {
      if (meshes.length === 1) {
        meshes[0].geometry.center();
        return;
      }

      var _getBoundingBox2 = this._getBoundingBox(meshes),
          min = _getBoundingBox2.min,
          max = _getBoundingBox2.max;

      this.center.addVectors(max, min).multiplyScalar(-0.5);
      this.offset.makeTranslation(this.center.x, this.center.y, this.center.z);

      for (var j = 0; j < meshes.length; j++) {
        var mesh = meshes[j];
        if (this.geometries.indexOf(mesh.geometry.uuid) !== -1) {
          continue;
        }
        mesh.geometry.applyMatrix(this.offset);
        this.geometries.push(mesh.geometry.uuid);
      }
    }
  }, {
    key: "_getBoundingBox",
    value: function _getBoundingBox(meshes) {
      for (var i = 0; i < meshes.length; ++i) {
        var mesh = meshes[i];
        if (!mesh.geometry.boundingBox) {
          mesh.geometry.computeBoundingBox();
        }
        var box = mesh.geometry.boundingBox;

        this.boundingBox.min.x = Math.min(box.min.x, box.min.x);
        this.boundingBox.min.y = Math.min(box.min.y, box.min.y);
        this.boundingBox.min.z = Math.min(box.min.z, box.min.z);

        this.boundingBox.max.x = Math.max(box.max.x, box.max.x);
        this.boundingBox.max.y = Math.max(box.max.y, box.max.y);
        this.boundingBox.max.z = Math.max(box.max.z, box.max.z);
      }
      return this.boundingBox;
    }
  }]);

  return ThreeToAmmo;
}();

module.exports = ThreeToAmmo;
