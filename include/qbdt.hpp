//////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano and the Qrack contributors 2017-2023. All rights reserved.
//
// QBinaryDecision tree is an alternative approach to quantum state representation, as
// opposed to state vector representation. This is a compressed form that can be
// operated directly on while compressed. Inspiration for the Qrack implementation was
// taken from JKQ DDSIM, maintained by the Institute for Integrated Circuits at the
// Johannes Kepler University Linz:
//
// https://github.com/iic-jku/ddsim
//
// Licensed under the GNU Lesser General Public License V3.
// See LICENSE.md in the project root or https://www.gnu.org/licenses/lgpl-3.0.en.html
// for details.

#pragma once

#include "mpsshard.hpp"
#include "qbdt_node.hpp"
#include "qengine.hpp"

#define QINTERFACE_TO_QALU(qReg) std::dynamic_pointer_cast<QAlu>(qReg)
#define QINTERFACE_TO_QPARITY(qReg) std::dynamic_pointer_cast<QParity>(qReg)

namespace Qrack {

class QBdt;
typedef std::shared_ptr<QBdt> QBdtPtr;

#if ENABLE_ALU
class QBdt : public QAlu, public QParity, public QInterface {
#else
class QBdt : public QParity, public QInterface {
#endif
protected:
    bitLenInt bdtStride;
    int64_t devID;
    QBdtNodeInterfacePtr root;
    bitCapInt bdtMaxQPower;
    std::vector<int64_t> deviceIDs;
    std::vector<QInterfaceEngine> engines;
    std::vector<MpsShardPtr> shards;

    void DumpBuffers() {
        for (size_t i = 0; i < shards.size(); ++i) {
            shards[i] = NULL;
        }
    }

    void FlushBuffer(bitLenInt t) {
        const MpsShardPtr shard = shards[t];
        if (shard) {
            shards[t] = NULL;
            ApplySingle(shard->gate, t);
        }
    }

    void FlushBuffers() {
        for (size_t i = 0U; i < shards.size(); ++i) {
            FlushBuffer(i);
        }
    }

    void FlushIfBlocked(bitLenInt target,
                        const std::vector<bitLenInt>& controls =
                                   std::vector<bitLenInt>()) {
      FlushIfBlocked(controls);
      FlushBuffer(target);
    }

    void FlushIfBlocked(const std::vector<bitLenInt>& controls) {
      for (const bitLenInt& control : controls) {
        const MpsShardPtr shard = shards[control];
        if (shard && !shard->IsPhase()) {
          shards[control] = NULL;
          ApplySingle(shard->gate, control);
        }
      }
    }

    void FlushNonPhaseBuffers() {
      for (size_t i = 0U; i < shards.size(); ++i) {
        const MpsShardPtr shard = shards[i];
        if (shard && !shard->IsPhase()) {
          shards[i] = NULL;
          ApplySingle(shard->gate, i);
        }
      }
    }

    QEnginePtr MakeQEngine(bitLenInt qbCount, bitCapInt perm = ZERO_BCI);

    template<typename Fn> void GetTraversal(Fn getLambda) {
      FlushBuffers();

      _par_for(maxQPower, [&] (const bitCapInt& i, const unsigned& cpu) {
               QBdtNodeInterfacePtr leaf = root;
               complex scale = leaf->scale;
               for (bitLenInt j = 0U; j < qubitCount; ++j) {
               leaf = leaf->branches[SelectBit(i, j)];
               if (!leaf) {
               break;
               }

               scale *= leaf->scale;
               }

               getLambda((bitCapIntOcl)i, scale);
               });
    }

    template<typename Fn> void SetTraversal(Fn setLambda) {
        DumpBuffers();
        root = std::make_shared<QBdtNode>();

#if ENABLE_QBDT_CPU_PARALLEL && ENABLE_PTHREAD
        if (true) {
          std::lock_guard<std::mutex> lock(root->mtx);
          root->Branch(qubitCount);
        }
#else
        root->Branch(qubitCount);
#endif

        _par_for(maxQPower, [&](const bitCapInt& i, const unsigned& cpu) {
                 QBdtNodeInterfacePtr leaf = root;
                 for (bitLenInt j = 0U; j < qubitCount; ++j) {
                 leaf = leaf->branches[SelectBit(i, j)];
                 }

                 setLambda((bitCapIntOcl)i, leaf);
                 });

        root->PopStateVector(qubitCount);
        root->Prune(qubitCount);
    }

    template<typename Fn> void ExecuteAsStateVector(Fn operation) {
      QInterfacePtr qReg = MakeQEngine(qubitCount);
      GetQuantumState(qReg);
      operation(qReg);
      SetQuantumState(qReg);
    }

    template<typename Fn> bitCapInt BitCapIntAsStateVector(Fn operation) {
      QInterfacePtr qReg = MakeQEngine(qubitCount);
      GetQuantumState(qReg);
      const bitCapInt toRet = operation(qReg);
      SetQuantumState(qReg);

      return toRet;
    }

    void par_for_qbdt(const bitCapInt& end, bitLenInt maxQubit, BdtFunc fn,
                      bool branch = true);

    void _par_for(const bitCapInt& end, ParallelFuncBdt fn);

    void DecomposeDispose(bitLenInt start, bitLenInt length, QBdtPtr dest);

    void ApplyControlledSingle(const complex* mtrx, std::vector<bitLenInt> controls,
                               bitLenInt target, bool isAnti);

    static size_t SelectBit(bitCapInt perm, bitLenInt bit) {
      return (size_t) bi_and_1(perm >> bit);
    }

    static bitCapInt RemovePower(bitCapInt perm, bitCapInt power) {
      bi_decrement(&power, 1U);
      return (perm & power) | ((perm >> 1U) & ~power);
    }

    void ApplySingle(const complex* mtrx, bitLenInt target);

    void Init();

    bitCapInt MAllOptionalCollapse(bool isCollapsing);

    bitCapInt SampleClone(const std::vector<bitCapInt>& qPowers) {
      const bitCapInt rawSample = MAllOptionalCollapse(false);
      bitCapInt sample = ZERO_BCI;
      for (size_t i = 0U; i < qPowers.size(); ++i) {
        if (bi_compare_0(rawSample & qPowers[i]) != 0) {
          bi_or_ip(&sample, pow2(i));
        }
      }

      return sample;
    }

public:
    QBdt(const std::vector<QInterfaceEngine>& eng, bitLenInt qBitCount,
         bitCapInt initState = ZERO_BCI, qrack_rand_gen_ptr rgp = nullptr,
         complex phaseFac = CMPLX_DEFAULT_ARG, bool doNorm = false,
         bool randomGlobalPhase = true, bool useHostMem = false,
         int64_t deviceId = -1, bool useHardwareRNG = true,
         bool useSparseStateVec = false, real1_f norm_thresh = REAL1_EPSILON,
         const std::vector<int64_t>& ignored = {},
         bitLenInt qubitThreshold = 0U,
         real1_f separation_thresh = FP_NORM_EPSILON_F);

    QBdt(bitLenInt qBitCount, bitCapInt initState = ZERO_BCI,
         qrack_rand_gen_ptr rgp = nullptr, complex phaseFac = CMPLX_DEFAULT_ARG,
         bool doNorm = false, bool randomGlobalPhase = true,
         bool useHostMem = false, int64_t deviceId = -1,
         bool useHardwareRNG = true, bool useSparseStateVec = false,
         real1_f norm_thresh = REAL1_EPSILON,
         const std::vector<int64_t>& devList = {},
         bitLenInt qubitThreshold = 0U,
         real1_f separation_thresh = FP_NORM_EPSILON_F)
      : QBdt({ QINTERFACE_OPTIMAL_BASE }, qBitCount, initState, rgp, phaseFac,
             doNorm, randomGlobalPhase, useHostMem, deviceId, useHardwareRNG,
             useSparseStateVec, norm_thresh, devList, qubitThreshold,
             separation_thresh) { }

    size_t CountBranches();

    bool isBinaryDecisionTree() { return true; };

    void SetDevice(int64_t dID) { devID = dID; }

    void UpdateRunningNorm(real1_f norm_thresh = REAL1_DEFAULT_ARG) { }

    void NormalizeState(real1_f nrm = REAL1_DEFAULT_ARG,
                        real1_f norm_thresh = REAL1_DEFAULT_ARG,
                        real1_f phaseArg = ZERO_R1_F) {
      root->Normalize(qubitCount);
    }

    real1_f SumSqrDiff(QInterfacePtr toCompare) {
      return SumSqrDiff(std::dynamic_pointer_cast<QBdt>(toCompare));
    }

    real1_f SumSqrDiff(QBdtPtr toCompare);

    void SetPermutation(bitCapInt initState,
                        const complex& phaseFac = CMPLX_DEFAULT_ARG);

    QInterfacePtr Clone();

    void GetQuantumState(complex* state) {
      GetTraversal([state](bitCapIntOcl i, complex scale) { state[i] = scale; });
    }

    void GetQuantumState(QInterfacePtr eng) {
      GetTraversal([eng](bitCapIntOcl i, complex scale) {
                   eng->SetAmplitude(i, scale);
                   });
    }

    void SetQuantumState(const complex* state) {
      SetTraversal([state](bitCapIntOcl i, QBdtNodeInterfacePtr leaf) {
                    leaf->scale = state[i];
                    });
    }

    void SetQuantumState(QInterfacePtr eng) {
      SetTraversal([eng](bitCapIntOcl i, QBdtNodeInterfacePtr leaf) {
                   leaf->scale = eng->GetAmplitude(i);
                   });
    }

    void GetProbs(real1* outputProbs) {
      GetTraversal([outputProbs](bitCapIntOcl i, complex scale) {
                   outputProbs[i] = norm(scale);
                   });
    }

    complex GetAmplitude(bitCapInt perm);

    void SetAmplitude(bitCapInt perm, const complex& amp) {
      ExecuteAsStateVector([&](QInterfacePtr eng) { eng->SetAmplitude(perm, amp); });
    }

    using QInterface::Compose;
    bitLenInt Compose(QBdtPtr toCopy, bitLenInt start);
    bitLenInt Compose(QInterfacePtr toCopy, bitLenInt start) {
      return Compose(std::dynamic_pointer_cast<QBdt>(toCopy), start);
    }

    void Decompose(bitLenInt start, QInterfacePtr dest) {
      QBdtPtr d = std::dynamic_pointer_cast<QBdt>(dest);
      DecomposeDispose(start, dest->GetQubitCount(), d);
    }

    QInterfacePtr Decompose(bitLenInt start, bitLenInt length);
    void Dispose(bitLenInt start, bitLenInt length) {
      DecomposeDispose(start, length, NULL);
    }

    void Dispose(bitLenInt start, bitLenInt length, bitCapInt disposedPerm) {
      ForceMReg(start, length, disposedPerm);
      DecomposeDispose(start, length, NULL);
    }

    using QInterface::Allocate;
    bitLenInt Allocate(bitLenInt start, bitLenInt length);

    real1_f Prob(bitLenInt qubitIndex);
    virtual real1_f ProbAll(bitCapInt fullRegister) override;

    bool ForceM(bitLenInt qubit, bool result, bool doForce = true, bool doApply = true);
    bitCapInt MAll() { return MAllOptionalCollapse(true); }

    void Mtrx(const complex* mtrx, bitLenInt target);
    void MCMtrx(const std::vector<bitLenInt>& controls, const complex* mtrx, bitLenInt target);
    void MACMtrx(const std::vector<bitLenInt>& controls, const complex* mtrx, bitLenInt target);
    void MCPhase(const std::vector<bitLenInt>& controls, complex topLeft, complex bottomRight, bitLenInt target);
    void MCInvert(const std::vector<bitLenInt>& controls, complex topRight, complex bottomLeft, bitLenInt target);

    void FSim(real1_f theta, real1_f phi, bitLenInt qubitIndex1, bitLenInt qubitIndex2);

    void Swap(bitLenInt q1, bitLenInt q2)
    {
      if (q2 < q1) {
        std::swap(q1, q2);
      }
      QInterface::Swap(q1, q2);
    }
    void ISwap(bitLenInt q1, bitLenInt q2)
    {
      if (q2 < q1) {
        std::swap(q1, q2);
      }
      QInterface::ISwap(q1, q2);
    }
    void IISwap(bitLenInt q1, bitLenInt q2)
    {
      if (q2 < q1) {
        std::swap(q1, q2);
      }
      QInterface::IISwap(q1, q2);
    }
    void SqrtSwap(bitLenInt q1, bitLenInt q2)
    {
      if (q2 < q1) {
        std::swap(q1, q2);
      }
      QInterface::SqrtSwap(q1, q2);
    }
    void ISqrtSwap(bitLenInt q1, bitLenInt q2)
    {
      if (q2 < q1) {
        std::swap(q1, q2);
      }
      QInterface::ISqrtSwap(q1, q2);
    }
    void CSwap(const std::vector<bitLenInt>& controls, bitLenInt q1, bitLenInt q2)
    {
      if (q2 < q1) {
        std::swap(q1, q2);
      }
      QInterface::CSwap(controls, q1, q2);
    }
    void CSqrtSwap(const std::vector<bitLenInt>& controls, bitLenInt q1, bitLenInt q2)
    {
      if (q2 < q1) {
        std::swap(q1, q2);
      }
      QInterface::CSqrtSwap(controls, q1, q2);
    }
    void CISqrtSwap(const std::vector<bitLenInt>& controls, bitLenInt q1, bitLenInt q2)
    {
      if (q2 < q1) {
        std::swap(q1, q2);
      }
      QInterface::CISqrtSwap(controls, q1, q2);
    }

    real1_f ProbParity(bitCapInt mask)
    {
      if (bi_compare_0(mask) == 0) {
        return ZERO_R1_F;
      }

      bitCapInt maskMin1 = mask;
      bi_decrement(&maskMin1, 1U);
      if (bi_compare_0(mask & maskMin1) == 0) {
        return Prob(log2(mask));
      }

      real1_f toRet;
      ExecuteAsStateVector([&](QInterfacePtr eng) { toRet = QINTERFACE_TO_QPARITY(eng)->ProbParity(mask); });

      return toRet;
    }
    void CUniformParityRZ(const std::vector<bitLenInt>& controls, bitCapInt mask, real1_f angle)
    {
      ExecuteAsStateVector(
        [&](QInterfacePtr eng) { QINTERFACE_TO_QPARITY(eng)->CUniformParityRZ(controls, mask, angle); });
    }
    bool ForceMParity(bitCapInt mask, bool result, bool doForce = true)
    {
      // If no bits in mask:
      if (bi_compare_0(mask) == 0) {
        return false;
      }

      // If only one bit in mask:
      bitCapInt maskMin1 = mask;
      bi_decrement(&maskMin1, 1U);
      if (bi_compare_0(mask & maskMin1) == 0) {
        return ForceM(log2(mask), result, doForce);
      }

      bool toRet;
      ExecuteAsStateVector(
        [&](QInterfacePtr eng) { toRet = QINTERFACE_TO_QPARITY(eng)->ForceMParity(mask, result, doForce); });

      return toRet;
    }

#if ENABLE_ALU
    using QInterface::M;
    bool M(bitLenInt q) { return QInterface::M(q); }
    using QInterface::X;
    void X(bitLenInt q) { QInterface::X(q); }
    void INC(bitCapInt toAdd, bitLenInt start, bitLenInt length) { QInterface::INC(toAdd, start, length); }
    void DEC(bitCapInt toSub, bitLenInt start, bitLenInt length) { QInterface::DEC(toSub, start, length); }
    void INCC(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt carryIndex)
    {
      QInterface::INCC(toAdd, start, length, carryIndex);
    }
    void DECC(bitCapInt toSub, bitLenInt start, bitLenInt length, bitLenInt carryIndex)
    {
      QInterface::DECC(toSub, start, length, carryIndex);
    }
    void INCS(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt overflowIndex)
    {
      QInterface::INCS(toAdd, start, length, overflowIndex);
    }
    void DECS(bitCapInt toSub, bitLenInt start, bitLenInt length, bitLenInt overflowIndex)
    {
      QInterface::DECS(toSub, start, length, overflowIndex);
    }
    void CINC(bitCapInt toAdd, bitLenInt inOutStart, bitLenInt length, const std::vector<bitLenInt>& controls)
    {
      QInterface::CINC(toAdd, inOutStart, length, controls);
    }
    void CDEC(bitCapInt toSub, bitLenInt inOutStart, bitLenInt length, const std::vector<bitLenInt>& controls)
    {
      QInterface::CDEC(toSub, inOutStart, length, controls);
    }
    void INCDECC(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt carryIndex)
    {
      QInterface::INCDECC(toAdd, start, length, carryIndex);
    }
    void MULModNOut(bitCapInt toMul, bitCapInt modN, bitLenInt inStart, bitLenInt outStart, bitLenInt length)
    {
      QInterface::MULModNOut(toMul, modN, inStart, outStart, length);
    }
    void IMULModNOut(bitCapInt toMul, bitCapInt modN, bitLenInt inStart, bitLenInt outStart, bitLenInt length)
    {
      QInterface::IMULModNOut(toMul, modN, inStart, outStart, length);
    }
    void CMULModNOut(bitCapInt toMul, bitCapInt modN, bitLenInt inStart, bitLenInt outStart, bitLenInt length,
                     const std::vector<bitLenInt>& controls)
    {
      QInterface::CMULModNOut(toMul, modN, inStart, outStart, length, controls);
    }
    void CIMULModNOut(bitCapInt toMul, bitCapInt modN, bitLenInt inStart, bitLenInt outStart, bitLenInt length,
                      const std::vector<bitLenInt>& controls)
    {
      QInterface::CIMULModNOut(toMul, modN, inStart, outStart, length, controls);
    }
    void PhaseFlipIfLess(bitCapInt greaterPerm, bitLenInt start, bitLenInt length)
    {
      ExecuteAsStateVector(
        [&](QInterfacePtr eng) { QINTERFACE_TO_QALU(eng)->PhaseFlipIfLess(greaterPerm, start, length); });
    }
    void CPhaseFlipIfLess(bitCapInt greaterPerm, bitLenInt start, bitLenInt length, bitLenInt flagIndex)
    {
      ExecuteAsStateVector([&](QInterfacePtr eng) {
                           QINTERFACE_TO_QALU(eng)->CPhaseFlipIfLess(greaterPerm, start, length, flagIndex);
                           });
    }
    void INCDECSC(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt overflowIndex, bitLenInt carryIndex)
    {
      ExecuteAsStateVector([&](QInterfacePtr eng) {
                           QINTERFACE_TO_QALU(eng)->INCDECSC(toAdd, start, length, overflowIndex, carryIndex);
                           });
    }
    void INCDECSC(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt carryIndex)
    {
      ExecuteAsStateVector(
        [&](QInterfacePtr eng) { QINTERFACE_TO_QALU(eng)->INCDECSC(toAdd, start, length, carryIndex); });
    }
#if ENABLE_BCD
    void INCBCD(bitCapInt toAdd, bitLenInt start, bitLenInt length)
    {
      ExecuteAsStateVector([&](QInterfacePtr eng) { QINTERFACE_TO_QALU(eng)->INCBCD(toAdd, start, length); });
    }
    void INCDECBCDC(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt carryIndex)
    {
      ExecuteAsStateVector(
        [&](QInterfacePtr eng) { QINTERFACE_TO_QALU(eng)->INCDECBCDC(toAdd, start, length, carryIndex); });
    }
#endif
    void MUL(bitCapInt toMul, bitLenInt inOutStart, bitLenInt carryStart, bitLenInt length)
    {
      ExecuteAsStateVector(
        [&](QInterfacePtr eng) { QINTERFACE_TO_QALU(eng)->MUL(toMul, inOutStart, carryStart, length); });
    }
    void DIV(bitCapInt toDiv, bitLenInt inOutStart, bitLenInt carryStart, bitLenInt length)
    {
      ExecuteAsStateVector(
        [&](QInterfacePtr eng) { QINTERFACE_TO_QALU(eng)->DIV(toDiv, inOutStart, carryStart, length); });
    }
    void POWModNOut(bitCapInt base, bitCapInt modN, bitLenInt inStart, bitLenInt outStart, bitLenInt length)
    {
      ExecuteAsStateVector(
        [&](QInterfacePtr eng) { QINTERFACE_TO_QALU(eng)->POWModNOut(base, modN, inStart, outStart, length); });
    }
    void CMUL(bitCapInt toMul, bitLenInt inOutStart, bitLenInt carryStart, bitLenInt length,
              const std::vector<bitLenInt>& controls)
    {
      ExecuteAsStateVector(
        [&](QInterfacePtr eng) { QINTERFACE_TO_QALU(eng)->CMUL(toMul, inOutStart, carryStart, length, controls); });
    }
    void CDIV(bitCapInt toDiv, bitLenInt inOutStart, bitLenInt carryStart, bitLenInt length,
              const std::vector<bitLenInt>& controls)
    {
      ExecuteAsStateVector(
        [&](QInterfacePtr eng) { QINTERFACE_TO_QALU(eng)->CDIV(toDiv, inOutStart, carryStart, length, controls); });
    }
    void CPOWModNOut(bitCapInt base, bitCapInt modN, bitLenInt inStart, bitLenInt outStart, bitLenInt length,
                     const std::vector<bitLenInt>& controls)
    {
      ExecuteAsStateVector([&](QInterfacePtr eng) {
                           QINTERFACE_TO_QALU(eng)->CPOWModNOut(base, modN, inStart, outStart, length, controls);
                           });
    }

    bitCapInt IndexedLDA(bitLenInt indexStart, bitLenInt indexLength,
                         bitLenInt valueStart, bitLenInt valueLength,
                         const void* values, bool resetValue = true) {
      return BitCapIntAsStateVector([&](QInterfacePtr eng) {
                                    return QINTERFACE_TO_QALU(eng)->IndexedLDA(
                                      indexStart, indexLength, valueStart,
                                      valueLength, values, resetValue); });
    }

    bitCapInt IndexedADC(bitLenInt indexStart, bitLenInt indexLength,
                         bitLenInt valueStart, bitLenInt valueLength,
                         bitLenInt carryIndex, const void* values) {
      return BitCapIntAsStateVector([&](QInterfacePtr eng) {
                                    return QINTERFACE_TO_QALU(eng)->IndexedADC(
                                      indexStart, indexLength, valueStart,
                                      valueLength, carryIndex, values); });
    }

    bitCapInt IndexedSBC(bitLenInt indexStart, bitLenInt indexLength,
                         bitLenInt valueStart, bitLenInt valueLength,
                         bitLenInt carryIndex, const void* values) {
      return BitCapIntAsStateVector([&](QInterfacePtr eng) {
                                    return QINTERFACE_TO_QALU(eng)->IndexedSBC(
                                      indexStart, indexLength, valueStart,
                                      valueLength, carryIndex, values); });
    }

    void Hash(bitLenInt start, bitLenInt length, const void* values) {
      ExecuteAsStateVector([&](QInterfacePtr eng) {
                           QINTERFACE_TO_QALU(eng)->Hash(start, length,
                                                         values); });
    }
#endif
};

} // namespace Qrack

