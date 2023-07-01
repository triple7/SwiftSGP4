// swift-tools-version: 5.8
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "SwiftSGP4",
    platforms: [
        .iOS("16"),
        .macOS("13"),
],
    products: [
        // Products define the executables and libraries a package produces, and make them visible to other packages.
        .library(
            name: "SwiftSGP4",
            targets: ["SGP4", "SwiftSGP4"]),
    ],
    dependencies: [
        .package(url: "https://github.com/triple7/SwiftCelesTrak", branch:"main"),
    ],
    targets: [
        .target(
            name: "SGP4",
            path: "Sources/SGP4"),
        .target(
            name: "SwiftSGP4",
            dependencies: ["SGP4", .product(name: "SwiftCelesTrak", package: "SwiftCelesTrak")])
    ]
)
