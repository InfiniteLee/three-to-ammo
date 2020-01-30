const assert = require("assert");
import * as THREE from "three";
import {
  iterateGeometries,
  TYPE,
  FIT,
  createCollisionShapes,
  createBoxShape,
  createCylinderShape,
  createCapsuleShape,
  createConeShape,
  createSphereShape,
  createHullShape,
  createHACDShapes,
  createVHACDShapes,
  createTriMeshShape,
  createHeightfieldTerrainShape
} from "../index";

import Ammo from "ammo.js/builds/ammo.wasm.js";
import AmmoWasm from "ammo.js/builds/ammo.wasm.wasm";
const AmmoModule = Ammo.bind(undefined, {
  locateFile(path) {
    if (path.endsWith(".wasm")) {
      return AmmoWasm;
    }
    return path;
  }
});

AmmoModule().then(Ammo => {
  describe("Three-to-Ammo API Tests", () => {
    context("helper methods", () => {
      it("iterateGeometries()", () => {
        "should extract all the vertices, matrices, and indexes from the mesh and submeshes",
          () => {
            const boxGeometry = new THREE.BoxBufferGeometry(1, 1, 1, 1, 1, 1);
            const boxMesh = new THREE.Mesh(boxGeometry);
            const vertices = [];
            const matrices = [];
            const indexes = [];
            iterateGeometries(boxMesh, {}, (vertexArray, matrixArray, indexArray) => {
              vertices.push(vertexArray);
              matrices.push(matrixArray);
              indexes.push(indexArray);
            });
            assert.equal(vertices.length, 1);
            assert.equal(matrices.length, 1);
            assert.equal(indexes.length, 1);
            assert.equal(vertices[0].length, 72);
            assert.equal(matrices[0].length, 16);
            assert.equal(indexes[0].length, 36);
          };
      });
    });

    describe("createCollisionShapes()", () => {
      const matrixWorld = new THREE.Matrix4();

      context("with options.fit === FIT.ALL", () => {
        const shapes = [
          TYPE.BOX,
          TYPE.CYLINDER,
          TYPE.CAPSULE,
          TYPE.CONE,
          TYPE.SPHERE,
          TYPE.HULL,
          TYPE.HACD,
          TYPE.VHACD,
          TYPE.MESH
        ];

        const vertices = [];
        const matrices = [];
        const indexes = [];

        const boxGeometry = new THREE.BoxBufferGeometry(1, 1, 1, 1, 1, 1);
        const boxMesh = new THREE.Mesh(boxGeometry);
        iterateGeometries(boxMesh, {}, (vertexArray, matrixArray, indexArray) => {
          vertices.push(vertexArray);
          matrices.push(matrixArray);
          indexes.push(indexArray);
        });

        for (let i = 0; i < shapes.length; i++) {
          const shape = shapes[i];
          context(`with options.type === TYPE.${shape.toUpperCase()}`, () => {
            it(`should create a ${shape} shape`, () => {
              const options = {
                type: shape,
                fit: FIT.ALL
              };
              const shapes = createCollisionShapes(vertices, matrices, indexes, matrixWorld.elements, options);
              assert.notEqual(shapes.length, 0);
              for (let j = 0; j < shapes.length; j++) {
                shapes[j].destroy();
              }
            });
          });
        }
      });
    });

    describe("createXXXShape()", () => {
      const matrixWorld = new THREE.Matrix4();

      context("with options.fit === FIT.MANUAL", () => {
        context("with options.type === TYPE.BOX", () => {
          it("should create a box shape", () => {
            const options = {
              fit: FIT.MANUAL,
              halfExtents: { x: 1, y: 1, z: 1 }
            };
            const shape = createBoxShape(null, null, matrixWorld.elements, options);
            assert.notEqual(shape, null);
            shape.destroy();
          });
        });

        context("with options.type === TYPE.CYLINDER", () => {
          it("should create a capsule shape", () => {
            const options = {
              fit: FIT.MANUAL,
              halfExtents: { x: 1, y: 1, z: 1 }
            };
            const shape = createCylinderShape(null, null, matrixWorld.elements, options);
            assert.notEqual(shape, null);
            shape.destroy();
          });
        });
        context("with options.type === TYPE.CAPSULE", () => {
          it("should create a capsule shape", () => {
            const options = {
              fit: FIT.MANUAL,
              halfExtents: { x: 1, y: 1, z: 1 }
            };
            const shape = createCapsuleShape(null, null, matrixWorld.elements, options);
            assert.notEqual(shape, null);
            shape.destroy();
          });
        });
        context("with options.type === TYPE.CONE", () => {
          it("should create a cone shape", () => {
            const options = {
              fit: FIT.MANUAL,
              halfExtents: { x: 1, y: 1, z: 1 }
            };
            const shape = createConeShape(null, null, matrixWorld.elements, options);
            assert.notEqual(shape, null);
            shape.destroy();
          });
        });
        context("with options.type === TYPE.SPHERE", () => {
          it("should create a sphere shape", () => {
            const options = {
              fit: FIT.MANUAL,
              sphereRadius: 1
            };
            const shape = createSphereShape(null, null, matrixWorld.elements, options);
            assert.notEqual(shape, null);
            shape.destroy();
          });
        });
      });

      context("with options.fit === FIT.ALL", () => {
        const vertices = [];
        const matrices = [];
        const indexes = [];
        before(() => {
          const boxGeometry = new THREE.BoxBufferGeometry(1, 1, 1, 1, 1, 1);
          const boxMesh = new THREE.Mesh(boxGeometry);
          iterateGeometries(boxMesh, {}, (vertexArray, matrixArray, indexArray) => {
            vertices.push(vertexArray);
            matrices.push(matrixArray);
            indexes.push(indexArray);
          });
        });
        context("with options.type === TYPE.BOX", () => {
          it("should create a box shape", () => {
            const shape = createBoxShape(vertices, matrices, matrixWorld.elements);
            assert.notEqual(shape, null);
            shape.destroy();
          });
        });
        context("with options.type === TYPE.CYLINDER", () => {
          it("should create a capsule shape", () => {
            const shape = createCylinderShape(vertices, matrices, matrixWorld.elements);
            assert.notEqual(shape, null);
            shape.destroy();
          });
        });
        context("with options.type === TYPE.CAPSULE", () => {
          it("should create a capsule shape", () => {
            const shape = createCapsuleShape(vertices, matrices, matrixWorld.elements);
            assert.notEqual(shape, null);
            shape.destroy();
          });
        });
        context("with options.type === TYPE.CONE", () => {
          it("should create a cone shape", () => {
            const shape = createConeShape(vertices, matrices, matrixWorld.elements);
            assert.notEqual(shape, null);
            shape.destroy();
          });
        });
        context("with options.type === TYPE.SPHERE", () => {
          it("should create a sphere shape", () => {
            const shape = createSphereShape(vertices, matrices, matrixWorld.elements);
            assert.notEqual(shape, null);
            shape.destroy();
          });
        });
        context("with options.type === TYPE.HULL", () => {
          it("should create a hull shape", () => {
            const shape = createHullShape(vertices, matrices, matrixWorld.elements);
            assert.notEqual(shape, null);
            shape.destroy();
          });
        });
        context("with options.type === TYPE.HACD", () => {
          it("should create at least one HACD shape", () => {
            const shapes = createHACDShapes(vertices, matrices, indexes, matrixWorld.elements);
            assert.notEqual(shapes.length, 0);
            for (let i = 0; i < shapes.length; i++) {
              shapes[i].destroy();
            }
          });
        });
        context("with options.type === TYPE.VHACD", () => {
          it("should create at least one VHACD shape", () => {
            const shapes = createVHACDShapes(vertices, matrices, indexes, matrixWorld.elements);
            assert.notEqual(shapes.length, 0);
            for (let i = 0; i < shapes.length; i++) {
              shapes[i].destroy();
            }
          });
        });
        context("with options.type === TYPE.MESH", () => {
          it("should create a mesh shape", () => {
            const shape = createTriMeshShape(vertices, matrices, indexes, matrixWorld.elements);
            assert.notEqual(shape, null);
            shape.destroy();
          });
        });
        context("with options.type === TYPE.HEIGHTFIELD", () => {
          it("should create a heightfield shape", () => {
            const options = {
              fit: FIT.MANUAL,
              heightfieldDistance: 1,
              heightfieldData: [
                [0, 0, 0, 0],
                [0, 1, 1, 0],
                [0, 1, 1, 0],
                [0, 0, 0, 0]
              ]
            };
            const shape = createHeightfieldTerrainShape(options);
            assert.notEqual(shape, null);
            shape.destroy();
          });
        });
      });
    });
  });
  run();
});
