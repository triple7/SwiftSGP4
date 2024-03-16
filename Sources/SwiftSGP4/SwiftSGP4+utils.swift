//
//  File.swift
//  
//
//  Created by Yuma decaux on 10/3/2024.
//

import Foundation
import SwiftCelesTrak

extension SwiftSGP4 {

    
    public func toggleActive(noradIds: [Int], active: Bool) {
        DispatchQueue.concurrentPerform(iterations: self.targetCount, execute:  { i in
            if noradIds.contains(self.targets[i].NORAD_CAT_ID) {
                self.targets[i].setActive(active: active)
            }
        })
    }
    
}
