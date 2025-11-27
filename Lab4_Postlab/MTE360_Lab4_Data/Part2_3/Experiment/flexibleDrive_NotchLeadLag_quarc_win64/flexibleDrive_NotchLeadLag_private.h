/*
 * flexibleDrive_NotchLeadLag_private.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "flexibleDrive_NotchLeadLag".
 *
 * Model version              : 1.75
 * Simulink Coder version : 9.6 (R2021b) 14-May-2021
 * C source code generated on : Sun Mar 19 18:48:03 2023
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_flexibleDrive_NotchLeadLag_private_h_
#define RTW_HEADER_flexibleDrive_NotchLeadLag_private_h_
#include "rtwtypes.h"
#include "multiword_types.h"
#include "zero_crossing_types.h"

/* A global buffer for storing error messages (defined in quanser_common library) */
EXTERN char _rt_error_message[512];
extern real_T rt_remd_snf(real_T u0, real_T u1);
extern real_T look1_binlxpw(real_T u0, const real_T bp0[], const real_T table[],
  uint32_T maxIndex);

/* private model entry point functions */
extern void flexibleDrive_NotchLeadLag_derivatives(void);

#endif                    /* RTW_HEADER_flexibleDrive_NotchLeadLag_private_h_ */
