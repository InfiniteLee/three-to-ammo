/* global Ammo,THREE */

const SHAPE_BOX = "box";
const SHAPE_CYLINDER = "cylinder";
const SHAPE_SPHERE = "sphere";
const SHAPE_CAPSULE = "capsule";
const SHAPE_CONE = "cone";
const SHAPE_HULL = "hull";
const SHAPE_MESH = "mesh";

class ThreeToAmmo {
  constructor() {
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

  createCollisionShape(sceneRoot, options) {
    const autoGenerateShape = options.autoGenerateShape || true;
    const shape = options.shape || SHAPE_HULL;
    const recenter = options.recenter || false;
    if (options.halfExtents) {
      this.halfExtents.set(options.halfExtents.x, options.halfExtents.y, options.halfExtents.z);
    }
    const cylinderAxis = options.cylinderAxis || "y";
    const sphereRadius = options.sphereRadius || 1;
    let collisionShape;
    let triMesh;
    let shapeHull;

    let meshes = this._getMeshes(sceneRoot);

    if (shape !== "mesh") {
      if (recenter) {
        this._recenter(sceneRoot, meshes);
      }
    }

    let vertices = this._getVertices(sceneRoot, meshes);
    boundingBox.setFromPoints(vertices);

    //TODO: Support convex hull decomposition, compound shapes, gimpact (dynamic trimesh)
    switch (shape) {
      case "box": {
        let { x, y, z } = this._getHalfExtents(boundingBox);
        let halfExtents = new Ammo.btVector3(x, y, z);
        collisionShape = new Ammo.btBoxShape(halfExtents);
        Ammo.destroy(halfExtents);
        break;
      }
      case "sphere": {
        let radius = 1;
        if (sphereRadius) {
          radius = sphereRadius;
        } else {
          let sphere = new THREE.Sphere();
          sphere.setFromPoints(vertices);
          if (isFinite(sphere.radius)) {
            radius = sphere.radius;
          }
        }
        collisionShape = new Ammo.btSphereShape(radius);
        break;
      }
      case "cylinder": {
        let { x, y, z } = this._getHalfExtents(boundingBox);
        let halfExtents = new Ammo.btVector3(x, y, z);
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
      case "capsule": {
        let { x, y, z } = this._getHalfExtents(boundingBox);
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
      case "cone": {
        let { x, y, z } = this._getHalfExtents(boundingBox);
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
      case "hull": {
        let scale = new Ammo.btVector3(sceneRoot.scale.x, sceneRoot.scale.y, sceneRoot.scale.z);
        let vec3 = new Ammo.btVector3();
        let originalHull = new Ammo.btConvexHullShape();
        originalHull.setMargin(data.margin);

        for (let i = 0; i < vertices.length; i++) {
          vec3.setValue(vertices[i].x, vertices[i].y, vertices[i].z);
          originalHull.addPoint(vec3, i == vertices.length - 1);
        }

        collisionShape = originalHull;
        if (originalHull.getNumVertices() >= 100) {
          //Bullet documentation says don't use convexHulls with 100 verts or more
          shapeHull = new Ammo.btShapeHull(originalHull);
          shapeHull.buildHull(data.margin);
          Ammo.destroy(originalHull);
          collisionShape = new Ammo.btConvexHullShape(
            Ammo.getPointer(shapeHull.getVertexPointer()),
            shapeHull.numVertices()
          );
        }
        collisionShape.setLocalScaling(scale);

        Ammo.destroy(scale);
        Ammo.destroy(vec3);
        break;
      }
      case "mesh": {
        if (data.type !== "static") {
          //TODO: support btTriangleMeshShape for dynamic trimeshes. (not performant)
          console.warn("non-static mesh colliders are not currently supported");
          break;
        }

        let a = new Ammo.btVector3();
        let b = new Ammo.btVector3();
        let c = new Ammo.btVector3();
        triMesh = new Ammo.btTriangleMesh(true, false);

        for (let j = 0; j < vertices.length; j += 3) {
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
      destroy: () => {
        if (shapeHull) Ammo.destroy(shapeHull);
        if (triMesh) Ammo.destroy(triMesh);
        if (collisionShape) Ammo.destroy(collisionShape);
      }
    };
  }

  _getMeshes(sceneRoot) {
    let meshes = [];
    sceneRoot.traverse(o => {
      if (o.type === "Mesh" && (!THREE.Sky || o.__proto__ != THREE.Sky.prototype)) {
        meshes.push(o);
      }
    });
    return meshes;
  }

  _getVertices(sceneRoot, meshes) {
    while (this.vertices.length > 0) {
      this.vertexPool.push(vertices.pop());
    }

    this.inverse.getInverse(sceneRoot.matrixWorld);

    for (let j = 0; j < meshes.length; j++) {
      let mesh = meshes[j];

      let geometry = mesh.geometry.index ? mesh.geometry.toNonIndexed() : mesh.geometry.clone();

      if (this.shape === SHAPE_MESH) {
        geometry.applyMatrix(mesh.matrixWorld);
      } else {
        this.matrix.multiplyMatrices(this.inverse, mesh.matrixWorld);
        geometry.applyMatrix(this.matrix);
      }

      if (geometry.isBufferGeometry) {
        let components = geometry.attributes.position.array;
        for (let i = 0; i < components.length; i += 3) {
          let x = components[i];
          let y = components[i + 1];
          let z = components[i + 2];

          if (this.vertexPool.length > 0) {
            this.vertices.push(this.vertexPool.pop().set(x, y, z));
          } else {
            this.vertices.push(new THREE.Vector3(x, y, z));
          }
        }
      } else {
        for (let i = 0; i < geometry.vertices.length; i++) {
          let vertex = geometry.vertices[i];
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

  _getHalfExtents(boundingBox) {
    if (this.autoGenerateShape) {
      let { min, max } = boundingBox;
      halfExtents.subVectors(max, min).multiplyScalar(0.5);
    }

    return {
      x: this.halfExtents.x,
      y: this.halfExtents.y,
      z: this.halfExtents.z
    };
  }

  _recenter(sceneRoot, meshes) {
    if (meshes.length === 1) {
      meshes[0].geometry.center();
      return;
    }

    let { min, max } = this._getBoundingBox(meshes);
    this.center.addVectors(max, min).multiplyScalar(-0.5);
    this.offset.makeTranslation(this.center.x, this.center.y, this.center.z);

    for (let j = 0; j < meshes.length; j++) {
      let mesh = meshes[j];
      if (this.geometries.indexOf(mesh.geometry.uuid) !== -1) {
        continue;
      }
      mesh.geometry.applyMatrix(this.offset);
      this.geometries.push(mesh.geometry.uuid);
    }
  }

  _getBoundingBox(meshes) {
    for (let i = 0; i < meshes.length; ++i) {
      let mesh = meshes[i];
      if (!mesh.geometry.boundingBox) {
        mesh.geometry.computeBoundingBox();
      }
      let box = mesh.geometry.boundingBox;

      this.boundingBox.min.x = Math.min(box.min.x, box.min.x);
      this.boundingBox.min.y = Math.min(box.min.y, box.min.y);
      this.boundingBox.min.z = Math.min(box.min.z, box.min.z);

      this.boundingBox.max.x = Math.max(box.max.x, box.max.x);
      this.boundingBox.max.y = Math.max(box.max.y, box.max.y);
      this.boundingBox.max.z = Math.max(box.max.z, box.max.z);
    }
    return this.boundingBox;
  }
}

module.exports = ThreeToAmmo;
