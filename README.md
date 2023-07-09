# SwiftSGP4

This is the swift wrapper for the SGP4 satellite orbit model. The SGP4 library is commonly used for predicting the positions of satellites based on orbital elements. The functions provide various utilities for initializing the satellite's orbit, propagating the orbit over time, converting coordinates, and performing calculations related to orbital elements and time.

the code was originally released in the 1980 and 1986
    spacetrack papers. a detailed discussion of the theory and history
    may be found in the [2006 aiaa paper](https://www.researchgate.net/publication/309393640_Orbit_Information_of_Predetermined_Accuracy_and_its_Sharing_in_the_SST_Context) by vallado, crawford, hujsak,
    and kelso.

                            companion code for
               fundamentals of astrodynamics and applications
                                    2013
                              by david vallado

     (w) 719-573-2600, email dvallado@agi.com, davallado@gmail.com

## features

This swift package module includes functions for initializing the satellite, propagating the orbit, converting between coordinate systems (e.g., Earth-centered Earth-fixed (ECEF) and geodetic), calculating satellite visibility (azimuth and elevation), and various mathematical operations for calculating between Date objects and JD epochs.

The single module dependency is [SwiftCelestrak](https://github.com/triple7/SwiftCelesTrak) which connects to the [celestrak satellite tracking service](https://celestrak.org) by Dr. T.S kelso.

The dependency is by the same author, Yuma Antoine Decaux, and is used for his own astronomical project, [AstreOS](https://astreos.space)

## Quick start

import SwiftCelesTrak
import SwiftSGP4

    func testSGP4propagation() {
        let GPGroups:[CelesTrakGroup] = [
            .active, // 7800
            ]
            let start = CACurrentMediaTime()
        self.celsTrak.getBatchGroupTargets(groups: GPGroups, returnFormat: .CSV, completion: { success in
            let end = CACurrentMediaTime()
            let count = self.celsTrak.targets.keys.count
            print("Batch of \(count) targets downloaded in \(end-start) seconds")
            
            // test propagating all downloaded targets
            let targets = Array(self.celsTrak.targets.values)
            let startSgp4 = CACurrentMediaTime()
            let sgp4 = SwiftSGP4(targets)
            let _ = sgp4.propagateOmms(1, 30)
            let endSgp4 = CACurrentMediaTime()
            print("propagated \(count) targets in \(endSgp4 - startSgp4) seconds")
        })
    }

This prints:
Batch of 8148 targets downloaded in 4.420734541665297 seconds
propagated 8148 targets in 0.03431679168716073 seconds
