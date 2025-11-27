/*
 * flexibleDrive_NotchLead.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "flexibleDrive_NotchLead".
 *
 * Model version              : 1.77
 * Simulink Coder version : 9.6 (R2021b) 14-May-2021
 * C source code generated on : Sun Mar 19 18:19:33 2023
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "flexibleDrive_NotchLead.h"
#include "flexibleDrive_NotchLead_private.h"
#include "flexibleDrive_NotchLead_dt.h"

/* Block signals (default storage) */
B_flexibleDrive_NotchLead_T flexibleDrive_NotchLead_B;

/* Continuous states */
X_flexibleDrive_NotchLead_T flexibleDrive_NotchLead_X;

/* Block states (default storage) */
DW_flexibleDrive_NotchLead_T flexibleDrive_NotchLead_DW;

/* Real-time model */
static RT_MODEL_flexibleDrive_NotchL_T flexibleDrive_NotchLead_M_;
RT_MODEL_flexibleDrive_NotchL_T *const flexibleDrive_NotchLead_M =
  &flexibleDrive_NotchLead_M_;
real_T look1_binlxpw(real_T u0, const real_T bp0[], const real_T table[],
                     uint32_T maxIndex)
{
  real_T frac;
  real_T yL_0d0;
  uint32_T bpIdx;
  uint32_T iLeft;
  uint32_T iRght;

  /* Column-major Lookup 1-D
     Search method: 'binary'
     Use previous index: 'off'
     Interpolation method: 'Linear point-slope'
     Extrapolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Linear'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u0 <= bp0[0U]) {
    iLeft = 0U;
    frac = (u0 - bp0[0U]) / (bp0[1U] - bp0[0U]);
  } else if (u0 < bp0[maxIndex]) {
    /* Binary Search */
    bpIdx = maxIndex >> 1U;
    iLeft = 0U;
    iRght = maxIndex;
    while (iRght - iLeft > 1U) {
      if (u0 < bp0[bpIdx]) {
        iRght = bpIdx;
      } else {
        iLeft = bpIdx;
      }

      bpIdx = (iRght + iLeft) >> 1U;
    }

    frac = (u0 - bp0[iLeft]) / (bp0[iLeft + 1U] - bp0[iLeft]);
  } else {
    iLeft = maxIndex - 1U;
    frac = (u0 - bp0[maxIndex - 1U]) / (bp0[maxIndex] - bp0[maxIndex - 1U]);
  }

  /* Column-major Interpolation 1-D
     Interpolation method: 'Linear point-slope'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  yL_0d0 = table[iLeft];
  return (table[iLeft + 1U] - yL_0d0) * frac + yL_0d0;
}

/*
 * This function updates continuous states using the ODE3 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
  static const real_T rt_ODE3_A[3] = {
    1.0/2.0, 3.0/4.0, 1.0
  };

  static const real_T rt_ODE3_B[3][3] = {
    { 1.0/2.0, 0.0, 0.0 },

    { 0.0, 3.0/4.0, 0.0 },

    { 2.0/9.0, 1.0/3.0, 4.0/9.0 }
  };

  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE3_IntgData *id = (ODE3_IntgData *)rtsiGetSolverData(si);
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T hB[3];
  int_T i;
  int_T nXc = 4;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  flexibleDrive_NotchLead_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  flexibleDrive_NotchLead_output();
  flexibleDrive_NotchLead_derivatives();

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  flexibleDrive_NotchLead_output();
  flexibleDrive_NotchLead_derivatives();

  /* tnew = t + hA(3);
     ynew = y + f*hB(:,3); */
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

real_T rt_remd_snf(real_T u0, real_T u1)
{
  real_T u1_0;
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1) || rtIsInf(u0)) {
    y = (rtNaN);
  } else if (rtIsInf(u1)) {
    y = u0;
  } else {
    if (u1 < 0.0) {
      u1_0 = ceil(u1);
    } else {
      u1_0 = floor(u1);
    }

    if ((u1 != 0.0) && (u1 != u1_0)) {
      u1_0 = fabs(u0 / u1);
      if (!(fabs(u1_0 - floor(u1_0 + 0.5)) > DBL_EPSILON * u1_0)) {
        y = 0.0 * u0;
      } else {
        y = fmod(u0, u1);
      }
    } else {
      y = fmod(u0, u1);
    }
  }

  return y;
}

/* Model output function */
void flexibleDrive_NotchLead_output(void)
{
  /* local block i/o variables */
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_NotchFilter;
  if (rtmIsMajorTimeStep(flexibleDrive_NotchLead_M)) {
    /* set solver stop time */
    if (!(flexibleDrive_NotchLead_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&flexibleDrive_NotchLead_M->solverInfo,
                            ((flexibleDrive_NotchLead_M->Timing.clockTickH0 + 1)
        * flexibleDrive_NotchLead_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&flexibleDrive_NotchLead_M->solverInfo,
                            ((flexibleDrive_NotchLead_M->Timing.clockTick0 + 1) *
        flexibleDrive_NotchLead_M->Timing.stepSize0 +
        flexibleDrive_NotchLead_M->Timing.clockTickH0 *
        flexibleDrive_NotchLead_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(flexibleDrive_NotchLead_M)) {
    flexibleDrive_NotchLead_M->Timing.t[0] = rtsiGetT
      (&flexibleDrive_NotchLead_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(flexibleDrive_NotchLead_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S1>/HIL Read Encoder Timebase' */

    /* S-Function Block: flexibleDrive_NotchLead/Plant/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder
        (flexibleDrive_NotchLead_DW.HILReadEncoderTimebase_Task, 1,
         &flexibleDrive_NotchLead_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        rtb_HILReadEncoderTimebase_o1 = 0;
        rtb_HILReadEncoderTimebase_o2 = 0;
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(flexibleDrive_NotchLead_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          flexibleDrive_NotchLead_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          flexibleDrive_NotchLead_DW.HILReadEncoderTimebase_Buffer[1];
      }
    }

    /* Gain: '<S1>/mm//cn' */
    flexibleDrive_NotchLead_B.mmcn = flexibleDrive_NotchLead_P.mmcn_Gain *
      rtb_HILReadEncoderTimebase_o1;

    /* S-Function (inverse_modulus_block): '<S1>/Inverse Modulus' */
    /* S-Function Block: flexibleDrive_NotchLead/Plant/Inverse Modulus (inverse_modulus_block) */
    {
      static const real_T sampling_period = 0.001;
      real_T half_range = flexibleDrive_NotchLead_P.InverseModulus_Modulus / 2.0;
      real_T du, dy;
      if (flexibleDrive_NotchLead_DW.InverseModulus_FirstSample) {
        flexibleDrive_NotchLead_DW.InverseModulus_FirstSample = false;
        flexibleDrive_NotchLead_DW.InverseModulus_PreviousInput =
          flexibleDrive_NotchLead_B.mmcn;
      }

      du = (real_T) flexibleDrive_NotchLead_B.mmcn -
        flexibleDrive_NotchLead_DW.InverseModulus_PreviousInput;
      if (du > half_range) {
        flexibleDrive_NotchLead_DW.InverseModulus_Revolutions =
          flexibleDrive_NotchLead_DW.InverseModulus_Revolutions - 1;
        dy = du - flexibleDrive_NotchLead_P.InverseModulus_Modulus;
      } else if (du < -half_range) {
        flexibleDrive_NotchLead_DW.InverseModulus_Revolutions =
          flexibleDrive_NotchLead_DW.InverseModulus_Revolutions + 1;
        dy = du + flexibleDrive_NotchLead_P.InverseModulus_Modulus;
      } else {
        dy = du;
      }

      flexibleDrive_NotchLead_B.InverseModulus = flexibleDrive_NotchLead_B.mmcn
        + flexibleDrive_NotchLead_DW.InverseModulus_Revolutions *
        flexibleDrive_NotchLead_P.InverseModulus_Modulus;
      flexibleDrive_NotchLead_DW.InverseModulus_PreviousInput =
        flexibleDrive_NotchLead_B.mmcn;
    }

    /* Gain: '<S1>/mm//cn - 2' */
    flexibleDrive_NotchLead_B.mmcn2 = flexibleDrive_NotchLead_P.mmcn2_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* S-Function (inverse_modulus_block): '<S1>/Inverse Modulus1' */
    /* S-Function Block: flexibleDrive_NotchLead/Plant/Inverse Modulus1 (inverse_modulus_block) */
    {
      static const real_T sampling_period = 0.001;
      real_T half_range = flexibleDrive_NotchLead_P.InverseModulus1_Modulus /
        2.0;
      real_T du, dy;
      if (flexibleDrive_NotchLead_DW.InverseModulus1_FirstSample) {
        flexibleDrive_NotchLead_DW.InverseModulus1_FirstSample = false;
        flexibleDrive_NotchLead_DW.InverseModulus1_PreviousInput =
          flexibleDrive_NotchLead_B.mmcn2;
      }

      du = (real_T) flexibleDrive_NotchLead_B.mmcn2 -
        flexibleDrive_NotchLead_DW.InverseModulus1_PreviousInput;
      if (du > half_range) {
        flexibleDrive_NotchLead_DW.InverseModulus1_Revolutions =
          flexibleDrive_NotchLead_DW.InverseModulus1_Revolutions - 1;
        dy = du - flexibleDrive_NotchLead_P.InverseModulus1_Modulus;
      } else if (du < -half_range) {
        flexibleDrive_NotchLead_DW.InverseModulus1_Revolutions =
          flexibleDrive_NotchLead_DW.InverseModulus1_Revolutions + 1;
        dy = du + flexibleDrive_NotchLead_P.InverseModulus1_Modulus;
      } else {
        dy = du;
      }

      flexibleDrive_NotchLead_B.InverseModulus1 =
        flexibleDrive_NotchLead_B.mmcn2 +
        flexibleDrive_NotchLead_DW.InverseModulus1_Revolutions *
        flexibleDrive_NotchLead_P.InverseModulus1_Modulus;
      flexibleDrive_NotchLead_DW.InverseModulus1_PreviousInput =
        flexibleDrive_NotchLead_B.mmcn2;
    }
  }

  /* SignalConversion: '<S2>/Output' incorporates:
   *  Clock: '<S2>/Clock'
   *  Constant: '<S2>/Constant'
   *  Lookup_n-D: '<S2>/Look-Up Table1'
   *  Math: '<S2>/Math Function'
   *  S-Function (sfun_tstart): '<S2>/startTime'
   *  Sum: '<S2>/Sum'
   */
  flexibleDrive_NotchLead_B.Output = look1_binlxpw(rt_remd_snf
    (flexibleDrive_NotchLead_M->Timing.t[0] - (0.0),
     flexibleDrive_NotchLead_P.Constant_Value),
    flexibleDrive_NotchLead_P.LookUpTable1_bp01Data, flexibleDrive_NotchLead_P.s,
    2000U);

  /* Gain: '<Root>/Gain' incorporates:
   *  Sum: '<Root>/Sum'
   */
  flexibleDrive_NotchLead_B.Gain = (flexibleDrive_NotchLead_B.Output -
    flexibleDrive_NotchLead_B.InverseModulus1) *
    flexibleDrive_NotchLead_P.Gain_Gain;

  /* TransferFcn: '<Root>/Lead Comp.1' */
  flexibleDrive_NotchLead_B.LeadComp1 = 0.0;
  flexibleDrive_NotchLead_B.LeadComp1 += flexibleDrive_NotchLead_P.LeadComp1_C *
    flexibleDrive_NotchLead_X.LeadComp1_CSTATE;
  flexibleDrive_NotchLead_B.LeadComp1 += flexibleDrive_NotchLead_P.LeadComp1_D *
    flexibleDrive_NotchLead_B.Gain;

  /* TransferFcn: '<Root>/Lead Comp.' */
  flexibleDrive_NotchLead_B.LeadComp = 0.0;
  flexibleDrive_NotchLead_B.LeadComp += flexibleDrive_NotchLead_P.LeadComp_C *
    flexibleDrive_NotchLead_X.LeadComp_CSTATE;
  flexibleDrive_NotchLead_B.LeadComp += flexibleDrive_NotchLead_P.LeadComp_D *
    flexibleDrive_NotchLead_B.LeadComp1;

  /* TransferFcn: '<Root>/Notch Filter' */
  rtb_NotchFilter = (flexibleDrive_NotchLead_P.NotchFilter_C[0] *
                     flexibleDrive_NotchLead_X.NotchFilter_CSTATE[0] +
                     flexibleDrive_NotchLead_P.NotchFilter_C[1] *
                     flexibleDrive_NotchLead_X.NotchFilter_CSTATE[1]) +
    flexibleDrive_NotchLead_P.NotchFilter_D * flexibleDrive_NotchLead_B.LeadComp;

  /* Saturate: '<Root>/Saturation' */
  if (rtb_NotchFilter > flexibleDrive_NotchLead_P.Saturation_UpperSat) {
    /* Saturate: '<Root>/Saturation' */
    flexibleDrive_NotchLead_B.Saturation =
      flexibleDrive_NotchLead_P.Saturation_UpperSat;
  } else if (rtb_NotchFilter < flexibleDrive_NotchLead_P.Saturation_LowerSat) {
    /* Saturate: '<Root>/Saturation' */
    flexibleDrive_NotchLead_B.Saturation =
      flexibleDrive_NotchLead_P.Saturation_LowerSat;
  } else {
    /* Saturate: '<Root>/Saturation' */
    flexibleDrive_NotchLead_B.Saturation = rtb_NotchFilter;
  }

  /* End of Saturate: '<Root>/Saturation' */

  /* Gain: '<S1>/Cable Gain' */
  flexibleDrive_NotchLead_B.CableGain = flexibleDrive_NotchLead_P.CableGain_Gain
    * flexibleDrive_NotchLead_B.Saturation;
  if (rtmIsMajorTimeStep(flexibleDrive_NotchLead_M)) {
    /* S-Function (hil_write_analog_block): '<S1>/HIL Write Analog' */

    /* S-Function Block: flexibleDrive_NotchLead/Plant/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      result = hil_write_analog(flexibleDrive_NotchLead_DW.HILInitialize_Card,
        &flexibleDrive_NotchLead_P.HILWriteAnalog_channels, 1,
        &flexibleDrive_NotchLead_B.CableGain);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(flexibleDrive_NotchLead_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void flexibleDrive_NotchLead_update(void)
{
  if (rtmIsMajorTimeStep(flexibleDrive_NotchLead_M)) {
    rt_ertODEUpdateContinuousStates(&flexibleDrive_NotchLead_M->solverInfo);
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++flexibleDrive_NotchLead_M->Timing.clockTick0)) {
    ++flexibleDrive_NotchLead_M->Timing.clockTickH0;
  }

  flexibleDrive_NotchLead_M->Timing.t[0] = rtsiGetSolverStopTime
    (&flexibleDrive_NotchLead_M->solverInfo);

  {
    /* Update absolute timer for sample time: [0.001s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++flexibleDrive_NotchLead_M->Timing.clockTick1)) {
      ++flexibleDrive_NotchLead_M->Timing.clockTickH1;
    }

    flexibleDrive_NotchLead_M->Timing.t[1] =
      flexibleDrive_NotchLead_M->Timing.clockTick1 *
      flexibleDrive_NotchLead_M->Timing.stepSize1 +
      flexibleDrive_NotchLead_M->Timing.clockTickH1 *
      flexibleDrive_NotchLead_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void flexibleDrive_NotchLead_derivatives(void)
{
  XDot_flexibleDrive_NotchLead_T *_rtXdot;
  _rtXdot = ((XDot_flexibleDrive_NotchLead_T *)
             flexibleDrive_NotchLead_M->derivs);

  /* Derivatives for TransferFcn: '<Root>/Lead Comp.1' */
  _rtXdot->LeadComp1_CSTATE = 0.0;
  _rtXdot->LeadComp1_CSTATE += flexibleDrive_NotchLead_P.LeadComp1_A *
    flexibleDrive_NotchLead_X.LeadComp1_CSTATE;
  _rtXdot->LeadComp1_CSTATE += flexibleDrive_NotchLead_B.Gain;

  /* Derivatives for TransferFcn: '<Root>/Lead Comp.' */
  _rtXdot->LeadComp_CSTATE = 0.0;
  _rtXdot->LeadComp_CSTATE += flexibleDrive_NotchLead_P.LeadComp_A *
    flexibleDrive_NotchLead_X.LeadComp_CSTATE;
  _rtXdot->LeadComp_CSTATE += flexibleDrive_NotchLead_B.LeadComp1;

  /* Derivatives for TransferFcn: '<Root>/Notch Filter' */
  _rtXdot->NotchFilter_CSTATE[0] = 0.0;
  _rtXdot->NotchFilter_CSTATE[0] += flexibleDrive_NotchLead_P.NotchFilter_A[0] *
    flexibleDrive_NotchLead_X.NotchFilter_CSTATE[0];
  _rtXdot->NotchFilter_CSTATE[1] = 0.0;
  _rtXdot->NotchFilter_CSTATE[0] += flexibleDrive_NotchLead_P.NotchFilter_A[1] *
    flexibleDrive_NotchLead_X.NotchFilter_CSTATE[1];
  _rtXdot->NotchFilter_CSTATE[1] +=
    flexibleDrive_NotchLead_X.NotchFilter_CSTATE[0];
  _rtXdot->NotchFilter_CSTATE[0] += flexibleDrive_NotchLead_B.LeadComp;
}

/* Model initialize function */
void flexibleDrive_NotchLead_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<S1>/HIL Initialize' */

  /* S-Function Block: flexibleDrive_NotchLead/Plant/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q2_usb", "0",
                      &flexibleDrive_NotchLead_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(flexibleDrive_NotchLead_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options
      (flexibleDrive_NotchLead_DW.HILInitialize_Card,
       "d0=digital;d1=digital;led=auto;update_rate=normal;decimation=1", 63);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(flexibleDrive_NotchLead_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(flexibleDrive_NotchLead_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(flexibleDrive_NotchLead_M, _rt_error_message);
      return;
    }

    if ((flexibleDrive_NotchLead_P.HILInitialize_AIPStart && !is_switching) ||
        (flexibleDrive_NotchLead_P.HILInitialize_AIPEnter && is_switching)) {
      flexibleDrive_NotchLead_DW.HILInitialize_AIMinimums[0] =
        (flexibleDrive_NotchLead_P.HILInitialize_AILow);
      flexibleDrive_NotchLead_DW.HILInitialize_AIMinimums[1] =
        (flexibleDrive_NotchLead_P.HILInitialize_AILow);
      flexibleDrive_NotchLead_DW.HILInitialize_AIMaximums[0] =
        flexibleDrive_NotchLead_P.HILInitialize_AIHigh;
      flexibleDrive_NotchLead_DW.HILInitialize_AIMaximums[1] =
        flexibleDrive_NotchLead_P.HILInitialize_AIHigh;
      result = hil_set_analog_input_ranges
        (flexibleDrive_NotchLead_DW.HILInitialize_Card,
         flexibleDrive_NotchLead_P.HILInitialize_AIChannels, 2U,
         &flexibleDrive_NotchLead_DW.HILInitialize_AIMinimums[0],
         &flexibleDrive_NotchLead_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(flexibleDrive_NotchLead_M, _rt_error_message);
        return;
      }
    }

    if ((flexibleDrive_NotchLead_P.HILInitialize_AOPStart && !is_switching) ||
        (flexibleDrive_NotchLead_P.HILInitialize_AOPEnter && is_switching)) {
      flexibleDrive_NotchLead_DW.HILInitialize_AOMinimums[0] =
        (flexibleDrive_NotchLead_P.HILInitialize_AOLow);
      flexibleDrive_NotchLead_DW.HILInitialize_AOMinimums[1] =
        (flexibleDrive_NotchLead_P.HILInitialize_AOLow);
      flexibleDrive_NotchLead_DW.HILInitialize_AOMaximums[0] =
        flexibleDrive_NotchLead_P.HILInitialize_AOHigh;
      flexibleDrive_NotchLead_DW.HILInitialize_AOMaximums[1] =
        flexibleDrive_NotchLead_P.HILInitialize_AOHigh;
      result = hil_set_analog_output_ranges
        (flexibleDrive_NotchLead_DW.HILInitialize_Card,
         flexibleDrive_NotchLead_P.HILInitialize_AOChannels, 2U,
         &flexibleDrive_NotchLead_DW.HILInitialize_AOMinimums[0],
         &flexibleDrive_NotchLead_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(flexibleDrive_NotchLead_M, _rt_error_message);
        return;
      }
    }

    if ((flexibleDrive_NotchLead_P.HILInitialize_AOStart && !is_switching) ||
        (flexibleDrive_NotchLead_P.HILInitialize_AOEnter && is_switching)) {
      flexibleDrive_NotchLead_DW.HILInitialize_AOVoltages[0] =
        flexibleDrive_NotchLead_P.HILInitialize_AOInitial;
      flexibleDrive_NotchLead_DW.HILInitialize_AOVoltages[1] =
        flexibleDrive_NotchLead_P.HILInitialize_AOInitial;
      result = hil_write_analog(flexibleDrive_NotchLead_DW.HILInitialize_Card,
        flexibleDrive_NotchLead_P.HILInitialize_AOChannels, 2U,
        &flexibleDrive_NotchLead_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(flexibleDrive_NotchLead_M, _rt_error_message);
        return;
      }
    }

    if (flexibleDrive_NotchLead_P.HILInitialize_AOReset) {
      flexibleDrive_NotchLead_DW.HILInitialize_AOVoltages[0] =
        flexibleDrive_NotchLead_P.HILInitialize_AOWatchdog;
      flexibleDrive_NotchLead_DW.HILInitialize_AOVoltages[1] =
        flexibleDrive_NotchLead_P.HILInitialize_AOWatchdog;
      result = hil_watchdog_set_analog_expiration_state
        (flexibleDrive_NotchLead_DW.HILInitialize_Card,
         flexibleDrive_NotchLead_P.HILInitialize_AOChannels, 2U,
         &flexibleDrive_NotchLead_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(flexibleDrive_NotchLead_M, _rt_error_message);
        return;
      }
    }

    if ((flexibleDrive_NotchLead_P.HILInitialize_EIPStart && !is_switching) ||
        (flexibleDrive_NotchLead_P.HILInitialize_EIPEnter && is_switching)) {
      flexibleDrive_NotchLead_DW.HILInitialize_QuadratureModes[0] =
        flexibleDrive_NotchLead_P.HILInitialize_EIQuadrature;
      flexibleDrive_NotchLead_DW.HILInitialize_QuadratureModes[1] =
        flexibleDrive_NotchLead_P.HILInitialize_EIQuadrature;
      result = hil_set_encoder_quadrature_mode
        (flexibleDrive_NotchLead_DW.HILInitialize_Card,
         flexibleDrive_NotchLead_P.HILInitialize_EIChannels, 2U,
         (t_encoder_quadrature_mode *)
         &flexibleDrive_NotchLead_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(flexibleDrive_NotchLead_M, _rt_error_message);
        return;
      }
    }

    if ((flexibleDrive_NotchLead_P.HILInitialize_EIStart && !is_switching) ||
        (flexibleDrive_NotchLead_P.HILInitialize_EIEnter && is_switching)) {
      flexibleDrive_NotchLead_DW.HILInitialize_InitialEICounts[0] =
        flexibleDrive_NotchLead_P.HILInitialize_EIInitial;
      flexibleDrive_NotchLead_DW.HILInitialize_InitialEICounts[1] =
        flexibleDrive_NotchLead_P.HILInitialize_EIInitial;
      result = hil_set_encoder_counts
        (flexibleDrive_NotchLead_DW.HILInitialize_Card,
         flexibleDrive_NotchLead_P.HILInitialize_EIChannels, 2U,
         &flexibleDrive_NotchLead_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(flexibleDrive_NotchLead_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S1>/HIL Read Encoder Timebase' */

  /* S-Function Block: flexibleDrive_NotchLead/Plant/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader
      (flexibleDrive_NotchLead_DW.HILInitialize_Card,
       flexibleDrive_NotchLead_P.HILReadEncoderTimebase_SamplesI,
       flexibleDrive_NotchLead_P.HILReadEncoderTimebase_Channels, 2,
       &flexibleDrive_NotchLead_DW.HILReadEncoderTimebase_Task);
    if (result >= 0) {
      result = hil_task_set_buffer_overflow_mode
        (flexibleDrive_NotchLead_DW.HILReadEncoderTimebase_Task,
         (t_buffer_overflow_mode)
         (flexibleDrive_NotchLead_P.HILReadEncoderTimebase_Overflow - 1));
    }

    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(flexibleDrive_NotchLead_M, _rt_error_message);
    }
  }

  /* Start for S-Function (inverse_modulus_block): '<S1>/Inverse Modulus' */

  /* S-Function Block: flexibleDrive_NotchLead/Plant/Inverse Modulus (inverse_modulus_block) */
  {
    flexibleDrive_NotchLead_DW.InverseModulus_FirstSample = true;
    flexibleDrive_NotchLead_DW.InverseModulus_Revolutions = 0;
  }

  /* Start for S-Function (inverse_modulus_block): '<S1>/Inverse Modulus1' */

  /* S-Function Block: flexibleDrive_NotchLead/Plant/Inverse Modulus1 (inverse_modulus_block) */
  {
    flexibleDrive_NotchLead_DW.InverseModulus1_FirstSample = true;
    flexibleDrive_NotchLead_DW.InverseModulus1_Revolutions = 0;
  }

  /* InitializeConditions for TransferFcn: '<Root>/Lead Comp.1' */
  flexibleDrive_NotchLead_X.LeadComp1_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<Root>/Lead Comp.' */
  flexibleDrive_NotchLead_X.LeadComp_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<Root>/Notch Filter' */
  flexibleDrive_NotchLead_X.NotchFilter_CSTATE[0] = 0.0;
  flexibleDrive_NotchLead_X.NotchFilter_CSTATE[1] = 0.0;
}

/* Model terminate function */
void flexibleDrive_NotchLead_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<S1>/HIL Initialize' */

  /* S-Function Block: flexibleDrive_NotchLead/Plant/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    hil_task_stop_all(flexibleDrive_NotchLead_DW.HILInitialize_Card);
    hil_monitor_stop_all(flexibleDrive_NotchLead_DW.HILInitialize_Card);
    is_switching = false;
    if ((flexibleDrive_NotchLead_P.HILInitialize_AOTerminate && !is_switching) ||
        (flexibleDrive_NotchLead_P.HILInitialize_AOExit && is_switching)) {
      flexibleDrive_NotchLead_DW.HILInitialize_AOVoltages[0] =
        flexibleDrive_NotchLead_P.HILInitialize_AOFinal;
      flexibleDrive_NotchLead_DW.HILInitialize_AOVoltages[1] =
        flexibleDrive_NotchLead_P.HILInitialize_AOFinal;
      num_final_analog_outputs = 2U;
    } else {
      num_final_analog_outputs = 0;
    }

    if (num_final_analog_outputs > 0) {
      result = hil_write_analog(flexibleDrive_NotchLead_DW.HILInitialize_Card,
        flexibleDrive_NotchLead_P.HILInitialize_AOChannels,
        num_final_analog_outputs,
        &flexibleDrive_NotchLead_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(flexibleDrive_NotchLead_M, _rt_error_message);
      }
    }

    hil_task_delete_all(flexibleDrive_NotchLead_DW.HILInitialize_Card);
    hil_monitor_delete_all(flexibleDrive_NotchLead_DW.HILInitialize_Card);
    hil_close(flexibleDrive_NotchLead_DW.HILInitialize_Card);
    flexibleDrive_NotchLead_DW.HILInitialize_Card = NULL;
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/

/* Solver interface called by GRT_Main */
#ifndef USE_GENERATED_SOLVER

void rt_ODECreateIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEDestroyIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEUpdateContinuousStates(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

#endif

void MdlOutputs(int_T tid)
{
  flexibleDrive_NotchLead_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  flexibleDrive_NotchLead_update();
  UNUSED_PARAMETER(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  flexibleDrive_NotchLead_initialize();
}

void MdlTerminate(void)
{
  flexibleDrive_NotchLead_terminate();
}

/* Registration function */
RT_MODEL_flexibleDrive_NotchL_T *flexibleDrive_NotchLead(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)flexibleDrive_NotchLead_M, 0,
                sizeof(RT_MODEL_flexibleDrive_NotchL_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&flexibleDrive_NotchLead_M->solverInfo,
                          &flexibleDrive_NotchLead_M->Timing.simTimeStep);
    rtsiSetTPtr(&flexibleDrive_NotchLead_M->solverInfo, &rtmGetTPtr
                (flexibleDrive_NotchLead_M));
    rtsiSetStepSizePtr(&flexibleDrive_NotchLead_M->solverInfo,
                       &flexibleDrive_NotchLead_M->Timing.stepSize0);
    rtsiSetdXPtr(&flexibleDrive_NotchLead_M->solverInfo,
                 &flexibleDrive_NotchLead_M->derivs);
    rtsiSetContStatesPtr(&flexibleDrive_NotchLead_M->solverInfo, (real_T **)
                         &flexibleDrive_NotchLead_M->contStates);
    rtsiSetNumContStatesPtr(&flexibleDrive_NotchLead_M->solverInfo,
      &flexibleDrive_NotchLead_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&flexibleDrive_NotchLead_M->solverInfo,
      &flexibleDrive_NotchLead_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&flexibleDrive_NotchLead_M->solverInfo,
      &flexibleDrive_NotchLead_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&flexibleDrive_NotchLead_M->solverInfo,
      &flexibleDrive_NotchLead_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&flexibleDrive_NotchLead_M->solverInfo,
                          (&rtmGetErrorStatus(flexibleDrive_NotchLead_M)));
    rtsiSetRTModelPtr(&flexibleDrive_NotchLead_M->solverInfo,
                      flexibleDrive_NotchLead_M);
  }

  rtsiSetSimTimeStep(&flexibleDrive_NotchLead_M->solverInfo, MAJOR_TIME_STEP);
  flexibleDrive_NotchLead_M->intgData.y = flexibleDrive_NotchLead_M->odeY;
  flexibleDrive_NotchLead_M->intgData.f[0] = flexibleDrive_NotchLead_M->odeF[0];
  flexibleDrive_NotchLead_M->intgData.f[1] = flexibleDrive_NotchLead_M->odeF[1];
  flexibleDrive_NotchLead_M->intgData.f[2] = flexibleDrive_NotchLead_M->odeF[2];
  flexibleDrive_NotchLead_M->contStates = ((real_T *) &flexibleDrive_NotchLead_X);
  rtsiSetSolverData(&flexibleDrive_NotchLead_M->solverInfo, (void *)
                    &flexibleDrive_NotchLead_M->intgData);
  rtsiSetSolverName(&flexibleDrive_NotchLead_M->solverInfo,"ode3");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = flexibleDrive_NotchLead_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;

    /* polyspace +2 MISRA2012:D4.1 [Justified:Low] "flexibleDrive_NotchLead_M points to
       static memory which is guaranteed to be non-NULL" */
    flexibleDrive_NotchLead_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    flexibleDrive_NotchLead_M->Timing.sampleTimes =
      (&flexibleDrive_NotchLead_M->Timing.sampleTimesArray[0]);
    flexibleDrive_NotchLead_M->Timing.offsetTimes =
      (&flexibleDrive_NotchLead_M->Timing.offsetTimesArray[0]);

    /* task periods */
    flexibleDrive_NotchLead_M->Timing.sampleTimes[0] = (0.0);
    flexibleDrive_NotchLead_M->Timing.sampleTimes[1] = (0.001);

    /* task offsets */
    flexibleDrive_NotchLead_M->Timing.offsetTimes[0] = (0.0);
    flexibleDrive_NotchLead_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(flexibleDrive_NotchLead_M,
             &flexibleDrive_NotchLead_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = flexibleDrive_NotchLead_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    flexibleDrive_NotchLead_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(flexibleDrive_NotchLead_M, 4.0);
  flexibleDrive_NotchLead_M->Timing.stepSize0 = 0.001;
  flexibleDrive_NotchLead_M->Timing.stepSize1 = 0.001;

  /* External mode info */
  flexibleDrive_NotchLead_M->Sizes.checksums[0] = (1467516827U);
  flexibleDrive_NotchLead_M->Sizes.checksums[1] = (3225263326U);
  flexibleDrive_NotchLead_M->Sizes.checksums[2] = (110936048U);
  flexibleDrive_NotchLead_M->Sizes.checksums[3] = (3196102782U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    flexibleDrive_NotchLead_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(flexibleDrive_NotchLead_M->extModeInfo,
      &flexibleDrive_NotchLead_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(flexibleDrive_NotchLead_M->extModeInfo,
                        flexibleDrive_NotchLead_M->Sizes.checksums);
    rteiSetTPtr(flexibleDrive_NotchLead_M->extModeInfo, rtmGetTPtr
                (flexibleDrive_NotchLead_M));
  }

  flexibleDrive_NotchLead_M->solverInfoPtr =
    (&flexibleDrive_NotchLead_M->solverInfo);
  flexibleDrive_NotchLead_M->Timing.stepSize = (0.001);
  rtsiSetFixedStepSize(&flexibleDrive_NotchLead_M->solverInfo, 0.001);
  rtsiSetSolverMode(&flexibleDrive_NotchLead_M->solverInfo,
                    SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  flexibleDrive_NotchLead_M->blockIO = ((void *) &flexibleDrive_NotchLead_B);
  (void) memset(((void *) &flexibleDrive_NotchLead_B), 0,
                sizeof(B_flexibleDrive_NotchLead_T));

  /* parameters */
  flexibleDrive_NotchLead_M->defaultParam = ((real_T *)
    &flexibleDrive_NotchLead_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &flexibleDrive_NotchLead_X;
    flexibleDrive_NotchLead_M->contStates = (x);
    (void) memset((void *)&flexibleDrive_NotchLead_X, 0,
                  sizeof(X_flexibleDrive_NotchLead_T));
  }

  /* states (dwork) */
  flexibleDrive_NotchLead_M->dwork = ((void *) &flexibleDrive_NotchLead_DW);
  (void) memset((void *)&flexibleDrive_NotchLead_DW, 0,
                sizeof(DW_flexibleDrive_NotchLead_T));

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    flexibleDrive_NotchLead_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 20;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* Initialize Sizes */
  flexibleDrive_NotchLead_M->Sizes.numContStates = (4);/* Number of continuous states */
  flexibleDrive_NotchLead_M->Sizes.numPeriodicContStates = (0);
                                      /* Number of periodic continuous states */
  flexibleDrive_NotchLead_M->Sizes.numY = (0);/* Number of model outputs */
  flexibleDrive_NotchLead_M->Sizes.numU = (0);/* Number of model inputs */
  flexibleDrive_NotchLead_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  flexibleDrive_NotchLead_M->Sizes.numSampTimes = (2);/* Number of sample times */
  flexibleDrive_NotchLead_M->Sizes.numBlocks = (25);/* Number of blocks */
  flexibleDrive_NotchLead_M->Sizes.numBlockIO = (10);/* Number of block outputs */
  flexibleDrive_NotchLead_M->Sizes.numBlockPrms = (4091);/* Sum of parameter "widths" */
  return flexibleDrive_NotchLead_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
