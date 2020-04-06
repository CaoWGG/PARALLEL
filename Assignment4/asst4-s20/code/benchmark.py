#!/usr/bin/python

import subprocess
import sys
import os
import os.path
import getopt
import math
import datetime
import random

import rutil

def usage(fname):
    
    ustring = "Usage: %s [-h][-g][-Q][-I][-b BENCHLIST] [-n NSTEP] [-p P1:P2:..:Pk] [-r RUNS] [-i ID] [-f OUTFILE]" % fname
    print ustring
    print "    -h            Print this message"
    print "    -g            Include mystery benchmarks for grading (Only available to graders)"
    print "    -Q            Quick mode: Don't compare with reference solution"
    print "    -I            Instrument simulation activities"
    print "    -b BENCHLIST  Specify which benchmark(s) to perform as substring of 'ABCDEF'"
    print "    -n NSTEP      Specify number of steps to run simulations"
    print "    -p P1:P2:..Pk Specify number of MPI processes as a colon-separated list"
    print "       If > 1, will run crun-mpi.  Else will run crun-seq"
    print "    -r RUNS       Set number of times each benchmark is run"
    print "    -i ID         Specify unique ID for distinguishing check files"
    print "    -f OUTFILE    Create output file recording measurements"
    print "         If file name contains field of form XX..X, will replace with ID having that many digits"
    sys.exit(0)

# General information
simProgram = "./crun-seq"
mpiSimProgram = "./crun-mpi"
refSimProgramDict = {'g': "./crun-soln-ghc", 'l': "./crun-soln-latedays", 'x' : ""}
defaultProcessCountsDict = { 'g': [8], 'l': [12], 'x' : [8] }

dataDirectory = "./data"

mpiFlagsDict = {'g': ["-map-by", "core", "-bind-to", "core"],
                'l': ["-bycore", "-bind-to-core"],
                'x': []}
outFile = None

doCheck = False
saveDirectory = "./check"

testFileName = ""
referenceFileName = ""

doInstrument = False


# How many times does each benchmark get run?
runCount = 3

# Grading parameters
pointsPerRun = 8
lowerThreshold = 0.60
upperThreshold = 0.95

# How many mismatched lines warrant detailed report
mismatchLimit = 5

# Graph/rat combinations: testId : (graphFile, ratFile, test name)
benchmarkDict = {
    'A': ('g-160x160-hlbrtA.gph', 'r-160x160-r40.rats', 'hlbrtA'),
    'B': ('g-160x160-hlbrtB.gph', 'r-160x160-r40.rats', 'hlbrtB'),
    'C': ('g-160x160-hlbrtC.gph', 'r-160x160-r40.rats', 'hlbrtC'),
    'D': ('g-160x160-hlbrtD.gph', 'r-160x160-r40.rats', 'hlbrtD'),
    'E': ('g-160x160-hlbrtE.gph', 'r-160x160-r40.rats', 'hlbrtE'),
    'F': ('g-160x160-hlbrtF.gph', 'r-160x160-r40.rats', 'hlbrtF'),
    }

graphWidth = 160
graphHeight = 160
loadFactor = 40
defaultSteps = 50
updateMode = "b"

defaultSeed = rutil.DEFAULTSEED

uniqueId = ""

def outmsg(s, noreturn = False):
    if len(s) > 0 and s[-1] != '\n' and not noreturn:
        s += "\n"
    sys.stdout.write(s)
    sys.stdout.flush()
    if outFile is not None:
        outFile.write(s)

def testName(testId, stepCount, seed, processCount):
    if testId in benchmarkDict:
        name = benchmarkDict[testId][-1]
    elif testId in extraBenchmarkDict:
        name = benchmarkDict[testId][-1]
    root = "%sx%.2d-step%.3d-seed%.3d" % (name, processCount, stepCount, seed)
    if uniqueId != "":
        root +=  ("-" + uniqueId)
    return root + ".txt"

def dataFile(fname):
    return dataDirectory + '/' + fname

def saveFileName(useRef, testId, stepCount, seed, processCount):
    return saveDirectory + "/" + ("ref" if useRef else "tst") + testName(testId, stepCount, seed, processCount)

def trim(s):
    while len(s) > 0 and s[-1] in '\r\n':
        s = s[:-1]
    return s

def checkOutputs(referenceFile, testFile, testName):
    if referenceFile == None or testFile == None:
        return True
    badLines = 0
    lineNumber = 0
    while True:
        rline = referenceFile.readline()
        tline = testFile.readline()
        lineNumber +=1
        if rline == "":
            if tline == "":
                break
            else:
                badLines += 1
                outmsg("Test %s.  Mismatch at line %d.  Reference simulation ended prematurely" % (testName, lineNumber))
                break
        elif tline == "":
            badLines += 1
            outmsg("Test %s.  Mismatch at line %d.  Simulation ended prematurely\n" % (testName, lineNumber))
            break
        rline = trim(rline)
        tline = trim(tline)
        if rline != tline:
            badLines += 1
            if badLines <= mismatchLimit:
                outmsg("Test %s.  Mismatch at line %d.  Expected result:'%s'.  Simulation result:'%s'\n" % (testName, lineNumber, rline, tline))
    referenceFile.close()
    testFile.close()
    if badLines > 0:
        outmsg("%d total mismatches.\n" % (badLines))
    return badLines == 0

def doRun(cmdList, simFileName):
    cmdLine = " ".join(cmdList)
    simFile = subprocess.PIPE
    if simFileName is not None:
        try:
            simFile = open(simFileName, 'w')
        except:
            print "Couldn't open output file '%s'" % fname
            return None
    tstart = datetime.datetime.now()
    try:
        outmsg("Running '%s'" % cmdLine)
        simProcess = subprocess.Popen(cmdList, stdout = simFile, stderr = subprocess.PIPE)
        simProcess.wait()
        if simFile != subprocess.PIPE:
            simFile.close()
        returnCode = simProcess.returncode
        # Echo any results printed by simulator on stderr onto stdout
        for line in simProcess.stderr:
            outmsg(line)
    except Exception as e:
        print "Execution of command '%s' failed. %s" % (cmdLine, e)
        if simFile != subprocess.PIPE:
            simFile.close()
        return None
    if returnCode == 0:
        delta = datetime.datetime.now() - tstart
        secs = delta.seconds + 24 * 3600 * delta.days + 1e-6 * delta.microseconds
        if simFile != subprocess.PIPE:
            simFile.close()
        return secs
    else:
        print "Execution of command '%s' gave return code %d" % (cmdLine, returnCode)
        if simFile != subprocess.PIPE:
            simFile.close()
        return None

def bestRun(cmdList, simFileName):
    sofar = 1e6
    for r in range(runCount):
        if runCount > 1:
            outmsg("Run #%d:" % (r+1), noreturn = True)
        secs = doRun(cmdList, simFileName)
        if secs is None:
            return None
        sofar = min(sofar, secs)
    return sofar

def runBenchmark(useRef, testId, stepCount, processCount, machine):
    global referenceFileName, testFileName
    nodes = graphWidth * graphHeight
    load = loadFactor
    gfname, rfname, tname = benchmarkDict[testId]
    graphFile = dataFile(gfname)
    ratFile = dataFile(rfname)
    params = [tname, str(stepCount)]
    results = params + [str(processCount)]
    preList = []
    prog = refSimProgramDict[machine] if useRef else simProgram if processCount == 1 else mpiSimProgram
    mpiFlags = mpiFlagsDict[machine]
    if processCount is None:
        processCount = processCountDict[machine]
    if processCount > 1:
        preList = ['mpirun', '-np', str(processCount)] + mpiFlags
    clist = ["-g", graphFile, "-r", ratFile, "-n", str(stepCount)]
    if doInstrument:
        clist += ["-I"]
    simFileName = None
    if not useRef:
        name = testName(testId, stepCount, defaultSeed, processCount)
        outmsg("+++++++++++++++++ Benchmark %s" % name)
    if doCheck:
        if not os.path.exists(saveDirectory):
            try:
                os.mkdir(saveDirectory)
            except Exception as e:
                outmsg("Couldn't create directory '%s' (%s)" % (saveDirectory, str(e)))
                simFile = subprocess.PIPE
        clist += ["-i", str(stepCount)]
        simFileName = saveFileName(useRef, testId, stepCount, defaultSeed, processCount)
        if useRef:
            referenceFileName = simFileName
        else:
            testFileName = simFileName
    else:
        clist += ["-q"]
    cmd = preList + [prog] + clist
    secs = bestRun(cmd, simFileName)
    if secs is None:
        return None
    else:
        rmoves = (nodes * load) * stepCount
        npm = 1e9 * secs/rmoves
        results.append("%.1f" % secs)
        results.append("%.1f" % npm)
        return results

def score(npm, rnpm):
    if npm == 0.0:
        return 0
    ratio = rnpm/npm
    nscore = 0.0
    if ratio >= upperThreshold:
        nscore = 1.0
    elif ratio >= lowerThreshold:
        nscore = (ratio-lowerThreshold)/(upperThreshold - lowerThreshold)
    # Round up to nearest 0.5
    return 0.5 * int(math.ceil(nscore * pointsPerRun * 2.0))

def formatTitle():
    ls = ["Name", "steps", "procs", "secs", "NPM"]
    if doCheck:
        ls += ["BNPM", "Ratio", "Pts"]
    return "\t".join(ls)

def sweep(testList, stepCount, processCount, machine):
    tcount = 0
    rcount = 0
    sum = 0.0
    refSum = 0.0
    resultList = []
    cresults = None
    totalPoints = 0
    for t in testList:
        ok = True
        results = runBenchmark(False, t, stepCount, processCount, machine)
        if results is not None and doCheck:
            cresults = runBenchmark(True, t, stepCount, processCount, machine)
            if referenceFileName != "" and testFileName != "":
                try:
                    rfile = open(referenceFileName, 'r')
                except:
                    rfile = None
                    print "Couldn't open reference simulation output file '%s'" % referenceFileName
                    ok = False
                try:
                    tfile = open(testFileName, 'r')
                except:
                    tfile = None
                    print "Couldn't open test simulation output file '%s'" % testFileName
                    ok = False
                if rfile is not None and tfile is not None:
                    ok = checkOutputs(rfile, tfile, t)
        if not ok:
            outmsg("TEST FAILED")
        if results is not None:
            tcount += 1
            npm = float(results[-1])
            sum += npm
            if cresults is not None:
                rcount += 1
                cnpm = float(cresults[-1])
                refSum += cnpm
                ratio = cnpm/npm if npm > 0 else 0
                points = score(npm, cnpm) if ok else 0.0
                totalPoints += points
                results += [cresults[-1], "%.2f" % ratio, "%.1f" % points]
            resultList.append(results)
    outmsg("+++++++++++++++++")
    outmsg(formatTitle())
    for r in resultList:
        outmsg("\t".join(r))
    if tcount > 0:
        avg = sum/tcount
        astring = "AVG:\t\t\t\t%.1f" % avg
        if refSum > 0:
            ravg = refSum/rcount
            astring += "\t%.1f" % ravg
        outmsg(astring)
        if doCheck:
            tstring = "TOTAL:\t\t\t\t\t\t\t%.1f" % totalPoints
            outmsg(tstring)

def generateFileName(template):
    global uniqueId
    myId = ""
    n = len(template)
    ls = []
    for i in range(n):
        c = template[i]
        if c == 'X':
            c = chr(random.randint(ord('0'), ord('9')))
        ls.append(c)
        myId += c
    if uniqueId == "":
        uniqueId = myId
    return "".join(ls) 

def run(name, args):
    global outFile, doCheck
    global uniqueId
    global runCount
    global doInstrument
    nstep = defaultSteps
    testList = None
    machine = 'x'
    try:
        host = os.environ['HOSTNAME']
    except:
        host = ''
    if host[:3] == 'ghc' or host[:4] == 'unix':
        machine = 'g'
        doCheck = True
    elif host[:8] == 'latedays' or host[:7] == 'compute':
        machine = 'l'
        doCheck = True
    else:
        outmsg("Warning: Host = '%s'. Can only get comparison results on GHC or Latedays machine" % host)
    processCounts = defaultProcessCountsDict[machine]

    optString = "hQIgb:n:p:r:i:f:"
    optlist, args = getopt.getopt(args, optString)
    for (opt, val) in optlist:
        if opt == '-h':
            usage(name)
        elif opt == '-Q':
            doCheck = False
        elif opt == '-I':
            doInstrument = True
        elif opt == '-g':
            try:
                import benchmark_extra
            except:
                print("Information about extra tests not available")
                return
            for k in benchmark_extra.extraBenchmarkDict.keys():
                benchmarkDict[k] = benchmark_extra.extraBenchmarkDict[k]
        elif opt == '-b':
            testList = list(val)
        elif opt == '-n':
            nstep = int(val)
        elif opt == '-r':
            runCount = int(val)
        elif opt == '-i':
            uniqueId = val
        elif opt == '-f':
            fname = generateFileName(val)
            try:
                outFile = open(fname, "w")
            except Exception as e:
                outFile = None
                outmsg("Couldn't open output file '%s'" % fname)
        elif opt == '-p':
            try:
                newProcessCounts = [int(s) for s in val.split(":")]
            except:
                print("Process counts must be given as colon-separated list")
                usage(name)
                return
            if min(newProcessCounts) < 1:
                print("Cannot have process count < 1")
                usage(name)
                return
            if max(newProcessCounts) > max(processCounts):
                print("Cannot have process count > %d" % max(processCounts))
                usage(name)
                return
            processCounts = newProcessCounts
        else:
            outmsg("Unknown option '%s'" % opt)
            usage(name)
    
    if testList is None:
        testList = sorted(benchmarkDict.keys())

    gstart = datetime.datetime.now()

    for p in processCounts:
        pstart = datetime.datetime.now()
        sweep(testList, nstep, p, machine)
        delta = datetime.datetime.now() - pstart
        secs = delta.seconds + 24 * 3600 * delta.days + 1e-6 * delta.microseconds
        print("Test time for %d processes = %.1f secs." % (p, secs))

    if len(processCounts) > 1:
        delta = datetime.datetime.now() - gstart
        secs = delta.seconds + 24 * 3600 * delta.days + 1e-6 * delta.microseconds
        print("Overall test time = %.1f secs." % (secs))

if __name__ == "__main__":
    run(sys.argv[0], sys.argv[1:])
