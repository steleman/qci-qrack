//////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano and the Qrack contributors 2017-2023. All rights reserved.
//
// This is a multithreaded, universal quantum register simulation, allowing
// (nonphysical) register cloning and direct measurement of probability and
// phase, to leverage what advantages classical emulation of qubits can have.
//
// Licensed under the GNU Lesser General Public License V3.
// See LICENSE.md in the project root or https://www.gnu.org/licenses/lgpl-3.0.en.html
// for details.

#include "qfactory.hpp"

namespace Qrack {

QBdtHybrid::QBdtHybrid(std::vector<QInterfaceEngine> eng, bitLenInt qBitCount, bitCapInt initState,
    qrack_rand_gen_ptr rgp, complex phaseFac, bool doNorm, bool randomGlobalPhase, bool useHostMem, int64_t deviceId,
    bool useHardwareRNG, bool useSparseStateVec, real1_f norm_thresh, std::vector<int64_t> devList,
    bitLenInt qubitThreshold, real1_f sep_thresh)
    : QInterface(qBitCount, rgp, doNorm, useHardwareRNG, randomGlobalPhase, doNorm ? norm_thresh : ZERO_R1_F)
    , useRDRAND(useHardwareRNG)
    , isSparse(useSparseStateVec)
    , useHostRam(useHostMem)
    , thresholdQubits(qubitThreshold)
    , separabilityThreshold(sep_thresh)
    , devID(deviceId)
    , engine(NULL)
    , phaseFactor(phaseFac)
    , deviceIDs(devList)
    , engines(eng)
{
    std::vector<QInterfaceEngine> e(engines);
    e.insert(e.begin(), QINTERFACE_BDT);
    qbdt = std::dynamic_pointer_cast<QBdt>(CreateQuantumInterface(e, qubitCount, initState, rand_generator, phaseFactor,
        doNormalize, randGlobalPhase, useHostRam, devID, useRDRAND, isSparse, (real1_f)amplitudeFloor, deviceIDs,
        thresholdQubits, separabilityThreshold));
}

QBdtHybrid::QBdtHybrid(QBdtPtr q, QEnginePtr e, std::vector<QInterfaceEngine> eng, bitLenInt qBitCount,
    bitCapInt initState, qrack_rand_gen_ptr rgp, complex phaseFac, bool doNorm, bool randomGlobalPhase, bool useHostMem,
    int64_t deviceId, bool useHardwareRNG, bool useSparseStateVec, real1_f norm_thresh, std::vector<int64_t> devList,
    bitLenInt qubitThreshold, real1_f sep_thresh)
    : QInterface(qBitCount, rgp, doNorm, useHardwareRNG, randomGlobalPhase, doNorm ? norm_thresh : ZERO_R1_F)
    , useRDRAND(useHardwareRNG)
    , isSparse(useSparseStateVec)
    , useHostRam(useHostMem)
    , thresholdQubits(qubitThreshold)
    , separabilityThreshold(sep_thresh)
    , devID(deviceId)
    , qbdt(q)
    , engine(e)
    , phaseFactor(phaseFac)
    , deviceIDs(devList)
    , engines(eng)
{
}

QInterfacePtr QBdtHybrid::MakeSimulator(bool isBdt, bitCapInt perm, complex phaseFac)
{
    std::vector<QInterfaceEngine> e(engines);
    e.insert(e.begin(), isBdt ? QINTERFACE_BDT : QINTERFACE_HYBRID);
    QInterfacePtr toRet =
        CreateQuantumInterface(e, 0U, perm, rand_generator, phaseFac, doNormalize, randGlobalPhase, useHostRam, devID,
            useRDRAND, isSparse, (real1_f)amplitudeFloor, deviceIDs, thresholdQubits, separabilityThreshold);
    std::dynamic_pointer_cast<QInterface>(toRet)->SetQubitCount(qubitCount);

    return toRet;
}

real1_f QBdtHybrid::ProbAll(bitCapInt fullRegister)
{
  const real1_f toRet = qbdt ?
                        qbdt->ProbAll(fullRegister) :
                        engine->ProbAll(fullRegister);

  if (toRet >= (ONE_R1_F - FP_NORM_EPSILON)) {
    SetPermutation(fullRegister);
  }

  return toRet;
}

} // namespace Qrack
