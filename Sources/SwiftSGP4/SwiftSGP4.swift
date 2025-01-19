import SwiftCelesTrak
import SGP4
import simd
import Foundation

public class SwiftSGP4 {
    private let pi:Double = 3.14159265358979323846
    private let jd1950:Double = 2433281.5
    private let rad:Double
    private let deg2rad:Double
    private let xpdotp:Double
    private let xpdotInv:Double
    private let xpdotInv2:Double
    private let minPDay:Double
    private let secPDay:Double
    
    public var targets:[CelesTrakTarget]
    private var satRecs:[elsetrec]
    public var coordinates:ContiguousArray<ContiguousArray<SIMD3<Double>>>
    // improved algorithm
    private let opsMode:CChar = "a".cString(using: .utf8)![0]
    // default second since last and fps
    private let secondsFromEpoch:Int = 1
    internal let fps:Int = 30

     public let targetCount:Int
    public let bufferCount:Int
    private let bufferOffset:Int
    internal var epochs = [Date]()
    internal var jdEpochs = [Double]()
    // last minute since value time propagation was calculated
        private var lastTimesSince = [Double]()
    internal var delta:Double = 0

    
    public init(_ targets: [CelesTrakTarget]) {
        
        self.targets = targets
//        let sattelliteSubsetTargets = [56768, 56771, 56773, 56775, 56776, 56782, 56783, 56786,
//                                       57885,57883,57882,57880,57881,57879,57878,57876,
//                                       57913, 57914, 57915, 57916, 57917, 57918, 57919, 57920, 57921]
//        self.targetCount = sattelliteSubsetTargets.count
        self.targetCount = targets.count
        self.rad = 180.0/self.pi
        self.deg2rad = pi / 180.0
        self.minPDay = 1440
        self.secPDay = 1440*60
        self.xpdotp = self.minPDay/(2.0*pi)
        self.xpdotInv = self.xpdotp*self.minPDay
        self.xpdotInv2 = self.xpdotp*self.minPDay*self.minPDay
        self.bufferCount = self.secondsFromEpoch*self.fps
        self.bufferOffset = self.bufferCount

        self.satRecs = [elsetrec]()
//        var target = targets.first!
//        print("jd \(jd)")
//        print(" jdfrac: \(jdFrac)")
//        print("bstar \(target.BSTAR)")
//        print("mean motion dot\(target.MEAN_MOTION_DOT)")
//        print("mean motion ddot\(target.MEAN_MOTION_DDOT)")
//        print("eccentricity\(target.ECCENTRICITY)")
//        print("argument of pericenter\(target.ARG_OF_PERICENTER)")
//        print("Inclination \(target.INCLINATION)")
//        print("mean anomaly \(target.MEAN_ANOMALY)")
//        print("mean motion\(target.MEAN_MOTION)")
//        print("rra of node\(target.RA_OF_ASC_NODE)")

        for target in targets {
//            if !sattelliteSubsetTargets.contains(target.NORAD_CAT_ID){
//                continue
//            }
            var satrec = elsetrec()
//            print(target.CLASSIFICATION_TYPE)
            satrec.noradID = Int32(target.NORAD_CAT_ID)
            satrec.elnum = target.ELEMENT_SET_NO
            satrec.revnum = target.REV_AT_EPOCH
            satrec.classification = target.CLASSIFICATION_TYPE.cString(using: .unicode)![0]
            satrec.ephtype = 0
            guard let epoch = dateString2Date(target.EPOCH) else {
                print("Skipping target with invalid EPOCH: \(target.EPOCH)")
                continue
            }
//            print("target \(target.NORAD_CAT_ID) has epoch at \(target.EPOCH)")
            let jdEpoch = timestampToJD(epoch)
            let currentDate = Date()
            let currentJd = timestampToJD(currentDate)
            
            let lastTSince = (currentJd - jdEpoch) * 1440.0
//            print("\(target.NORAD_CAT_ID) lastTSince: \(lastTSince)")
            
            _ = sgp4init(wgs72, opsMode, &genSatNum
                         , jdEpoch - jd1950, target.BSTAR, target.MEAN_MOTION_DOT/xpdotInv, target.MEAN_MOTION_DDOT/xpdotInv2, target.ECCENTRICITY, target.ARG_OF_PERICENTER*deg2rad, target.INCLINATION*deg2rad, target.MEAN_ANOMALY*deg2rad,
                         target.MEAN_MOTION/xpdotp, target.RA_OF_ASC_NODE*deg2rad, &satrec)

            
//            // sgp4fix allow additional parameters in the struct
//            satrec->ndot = xndot;
//            satrec->nddot = xnddot;
//            satrec->ecco = xecco;
//            satrec->argpo = xargpo;
//            satrec->inclo = xinclo;
//            satrec->mo = xmo;
//            // sgp4fix rename variables to clarify which mean motion is intended
//            satrec->no_kozai = xno_kozai;
//            satrec->nodeo = xnodeo;            
            satRecs.append(satrec)
            epochs.append(epoch)
            jdEpochs.append(jdEpoch)
            lastTimesSince.append(lastTSince)
            
        }

        let targetFrames = ContiguousArray<SIMD3<Double>>(repeating: zeroSimd, count: self.bufferCount + self.bufferOffset)
        self.coordinates = ContiguousArray<ContiguousArray<SIMD3<Double>>>(repeating: targetFrames, count: targetCount)
    }

    fileprivate var lastAppTSince:Double = 0
    public func propagateOmmsSingle(){
        // calculate the JD timestamp right now to pass into the ITRF function
        let currentDate = Date()
        let currentJd = timestampToJD(currentDate)

        DispatchQueue.concurrentPerform(iterations: self.targetCount, execute:  { satrecIndex in
            var satrec = satRecs[satrecIndex]
            // Calculate the target states from epoch to secondsFromEpoch
            let lastSince:Double = (currentJd - jdEpochs[satrecIndex]) * 1440.0
            var ro = [Double](repeating: 0, count: 3)
            var vo = [Double](repeating: 0, count: 3)
            sgp4(&satrec, lastSince, &ro, &vo)
            self.coordinates[satrecIndex][0] = SIMD3<Double>(ro)
        })
    }
    
    public func propagateOmmsByTargetJds(_ targetJds: [Double], _ noradID: Int) -> [[Double]] {
        // Initialize an array to store the return positions for each target JD
        var returnPositions: [[Double]] = Array(repeating: [], count: targetJds.count)
        for satrecIndex in 0..<self.targetCount {
            var satrec = satRecs[satrecIndex]
            // Check if the current satellite matches the specified NORAD ID
            if satrec.noradID == noradID {
                for (i, targetJd) in targetJds.enumerated() {
                    // Calculate the elapsed time since epoch in minutes
                    let elapsedTimeSinceEpoch: Double = (targetJd - jdEpochs[satrecIndex]) * 1440.0
                    // Prepare arrays for position (ro) and velocity (vo)
                    var ro = [Double](repeating: 0, count: 3)
                    var vo = [Double](repeating: 0, count: 3)
                    // Call sgp4 to compute the position and velocity
                    sgp4(&satrec, elapsedTimeSinceEpoch, &ro, &vo)
                    // Store the calculated position (ro) in the result
                    returnPositions[i] = ro
                }
                // Break after finding the matching NORAD ID
                break
            }
        }
        
        return returnPositions
    }

    public func propagateOmmsMinutes(_ minutes: Int = 60, _ noradID: Int) -> [[Double]] {
        // Get the JD timestamp of the current timestamp then go through the next minutes to create a set of points that the satellite will take
        let currentDate = Date()
        let currentJd = timestampToJD(currentDate)
        var returnPositions: [[Double]] = Array(repeating: [], count: minutes)
        for satrecIndex in 0..<self.targetCount{
            var satrec = satRecs[satrecIndex]
            if satrec.noradID == noradID{
                // Calculate the orbit positions for n minutes for this satellite
                for i in 0..<minutes{
                    // Calculate the target states from epoch to secondsFromEpoch
                    let lastSince:Double = (currentJd - jdEpochs[satrecIndex]) * 1440.0 + Double(i)
                    var ro = [Double](repeating: 0, count: 3)
                    var vo = [Double](repeating: 0, count: 3)
                    sgp4(&satrec, lastSince, &ro, &vo)
                    returnPositions[i] = ro
                }
            }
        }
        
        return returnPositions
    }

    public func propagateOmmsByDateTimestamp(dates: [Date]) -> [[SIMD3<Double>]] {
        var output = [[SIMD3<Double>]](repeating: [SIMD3<Double>](repeating: SIMD3<Double>(x: 0, y: 0, z: 0), count: dates.count), count: self.targetCount)
        // calculate the JD timestamp for each date to calculate the ITRF
        for (i, date) in dates.enumerated() {
            let currentJd = timestampToJD(date)
            
            DispatchQueue.concurrentPerform(iterations: self.targetCount, execute:  { satrecIndex in
                var satrec = satRecs[satrecIndex]
                // Calculate the target states from epoch to secondsFromEpoch
                let lastSince:Double = (currentJd - jdEpochs[satrecIndex]) * 1440.0
                var ro = [Double](repeating: 0, count: 3)
                var vo = [Double](repeating: 0, count: 3)
                sgp4(&satrec, lastSince, &ro, &vo)
                output[satrecIndex][i] = SIMD3<Double>(ro)
            })
        }
            return output
    }

    public func propagateOmms( _ minDelta: Double = 1/60 /* seconds */) {
        
        // time dimension parameters
        // We are propagating from
        // epoch to secondsFromEpoch by frames per second
        // such that:
        // delta = (1/(60*fps)
        // and count = seconds*fps
        delta = 1/Double(secondsFromEpoch*60*fps)

        // calculate the JD timestamp right now to pass into the ITRF function
        let currentDate = Date()
        let currentJd = timestampToJD(currentDate)

        DispatchQueue.concurrentPerform(iterations: self.targetCount, execute:  { i in
            computeITRF(i, jdEpochs[i], delta, currentJd)
        })
        // Double buffer to cycle around
        self.currentBufferOffset = self.currentBufferOffset + self.bufferOffset
        if self.currentBufferOffset == self.bufferCount*2 {
            self.currentBufferOffset = 0
        }
        lastAppTSince += Double(fps)*delta
    }
    
    private let zeroSimd = SIMD3<Double>([0, 0, 0])
    // Using generic number as we don't need satNum for propagation
    private var genSatNum = (Int8(77), Int8(77), Int8(77), Int8(77), Int8(77), Int8(77))
    // using index to circle around the frames buffer
    private var currentBufferOffset:Int = 0
    public func computeITRF(_ satrecIndex: Int, _ epoch: Double, _ delta: Double, _ currentJd: Double) {
        // struct to pass to sgp4 function
        var satrec = satRecs[satrecIndex]
        // Calculate the target states from epoch to secondsFromEpoch
        var lastSince:Double = 0
        lastSince = (currentJd - epoch) * 1440.0
        
        DispatchQueue.concurrentPerform(iterations: self.bufferCount, execute:  { i in
            // orbital set
            var ro = [Double](repeating: 0, count: 3)
            var vo = [Double](repeating: 0, count: 3)

            lastSince = lastSince + Double(i+1)*delta/15.0
            
            sgp4(&satrec, lastSince, &ro, &vo)
//            // transform from TEME to GTRF
//            var RGtrf = [Double](repeating: 0, count: 3)
//            let gmst = gstime(jdut1: epoch)
//            let gmstCos = cos(gmst)
//            let gmstSin = sin(gmst)
////
//            teme2ecefOptimised(&ro, epoch, gmstCos, gmstSin, &RGtrf)
//            self.coordinates[satrecIndex][i + currentBufferOffset] = SIMD3<Double>(RGtrf)
            self.coordinates[satrecIndex][i + currentBufferOffset] = SIMD3<Double>(ro)
        })
    }
    
}
