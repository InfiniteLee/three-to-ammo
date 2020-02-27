const path = require("path");
const webpack = require("webpack");

module.exports = (env, argv) => ({
  node: {
    fs: "empty"
  },
  entry: { test: "./test/test.js" },
  output: {
    path: path.resolve(__dirname, "test"),
    publicPath: "/test/",
    filename: "[name].js"
  },
  plugins: [new webpack.ProvidePlugin({ THREE: "three" })],
  devtool: "inline-source-map",
  devServer: {
    host: "0.0.0.0"
  },
  module: {
    rules: [
      {
        test: /\.(wasm)$/,
        type: "javascript/auto",
        use: {
          loader: "file-loader",
          options: {
            outputPath: "dist",
            name: "[name]-[hash].[ext]"
          }
        }
      }
    ]
  }
});
