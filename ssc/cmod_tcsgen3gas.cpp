/*******************************************************************************************************
*  Copyright 2017 Alliance for Sustainable Energy, LLC
*
*  NOTICE: This software was developed at least in part by Alliance for Sustainable Energy, LLC
*  ("Alliance") under Contract No. DE-AC36-08GO28308 with the U.S. Department of Energy and the U.S.
*  The Government retains for itself and others acting on its behalf a nonexclusive, paid-up,
*  irrevocable worldwide license in the software to reproduce, prepare derivative works, distribute
*  copies to the public, perform publicly and display publicly, and to permit others to do so.
*
*  Redistribution and use in source and binary forms, with or without modification, are permitted
*  provided that the following conditions are met:
*
*  1. Redistributions of source code must retain the above copyright notice, the above government
*  rights notice, this list of conditions and the following disclaimer.
*
*  2. Redistributions in binary form must reproduce the above copyright notice, the above government
*  rights notice, this list of conditions and the following disclaimer in the documentation and/or
*  other materials provided with the distribution.
*
*  3. The entire corresponding source code of any redistribution, with or without modification, by a
*  research entity, including but not limited to any contracting manager/operator of a United States
*  National Laboratory, any institution of higher learning, and any non-profit organization, must be
*  made publicly available under this license for as long as the redistribution is made available by
*  the research entity.
*
*  4. Redistribution of this software, without modification, must refer to the software by the same
*  designation. Redistribution of a modified version of this software (i) may not refer to the modified
*  version by the same designation, or by any confusingly similar designation, and (ii) must refer to
*  the underlying software originally provided by Alliance as �System Advisor Model� or �SAM�. Except
*  to comply with the foregoing, the terms �System Advisor Model�, �SAM�, or any confusingly similar
*  designation may not be used to refer to any modified version of this software or any modified
*  version of the underlying software originally provided by Alliance without the prior written consent
*  of Alliance.
*
*  5. The name of the copyright holder, contributors, the United States Government, the United States
*  Department of Energy, or any of their employees may not be used to endorse or promote products
*  derived from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
*  IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
*  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER,
*  CONTRIBUTORS, UNITED STATES GOVERNMENT OR UNITED STATES DEPARTMENT OF ENERGY, NOR ANY OF THEIR
*  EMPLOYEES, BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
*  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
*  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
*  IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
*  THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*******************************************************************************************************/

#include "core.h"

// for adjustment factors
#include "common.h"

// solarpilot header files
#include "AutoPilot_API.h"
#include "SolarField.h"
#include "IOUtil.h"
#include "csp_common.h"

// Can probably delete these headers later...
#include "csp_solver_util.h"
#include "csp_solver_core.h"
#include "csp_solver_pt_sf_perf_interp.h"
#include "csp_solver_gas_receiver.h"
#include "csp_solver_gas_collector_receiver.h"
#include "csp_solver_pc_Rankine_indirect_224.h"		// include for user-defined cycle
#include "csp_solver_pc_sco2.h"
#include "csp_solver_two_tank_tes.h"
#include "csp_solver_tou_block_schedules.h"

#include "csp_system_costs.h"

static var_info _cm_vtab_gen3gas[] = {
    // VARTYPE       DATATYPE    NAME                                  LABEL                                                                                                                                      UNITS           META                                 GROUP                                       REQUIRED_IF                                                         CONSTRAINTS      UI_HINTS
    { SSC_INPUT,     SSC_STRING, "solar_resource_file",                "Local weather file path",                                                                                                                 "",             "",                                  "Location and Resource",                    "?",                                                                "LOCAL_FILE",    ""},
    { SSC_INPUT,     SSC_TABLE,  "solar_resource_data",                "Weather resource data in memory",                                                                                                         "",             "",                                  "Location and Resource",                    "?",                                                                "",              ""},

    { SSC_INPUT,     SSC_NUMBER, "ppa_multiplier_model",               "PPA multiplier model",                                                                                                                    "0/1",          "0=diurnal,1=timestep",              "Time of Delivery Factors",                 "?=0",                                                              "INTEGER,MIN=0", ""},
    { SSC_INPUT,     SSC_ARRAY,  "dispatch_factors_ts",                "Dispatch payment factor array",                                                                                                           "",             "",                                  "Time of Delivery Factors",                 "ppa_multiplier_model=1",                                           "",              ""},

    { SSC_INPUT,     SSC_NUMBER, "field_model_type",                   "3=user performance maps vs solar position",           "",             "",                                  "Heliostat Field",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "gross_net_conversion_factor",        "Estimated gross to net conversion factor",                                                                                                "",             "",                                  "System Design",                            "*",                                                                "",              ""},

    { SSC_INPUT,     SSC_NUMBER, "dni_des",                            "Design-point DNI",                                                                                                                        "W/m2",         "",                                  "System Design",                            "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "p_start",                            "Heliostat startup energy",                                                                                                                "kWe-hr",       "",                                  "Heliostat Field",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "p_track",                            "Heliostat tracking energy",                                                                                                               "kWe",          "",                                  "Heliostat Field",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "hel_stow_deploy",                    "Stow/deploy elevation angle",                                                                                                             "deg",          "",                                  "Heliostat Field",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "v_wind_max",                         "Heliostat max wind velocity",                                                                                                             "m/s",          "",                                  "Heliostat Field",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "interp_nug",                         "Interpolation nugget",                                                                                                                    "-",            "",                                  "Heliostat Field",                          "?=0",                                                              "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "interp_beta",                        "Interpolation beta coef.",                                                                                                                "-",            "",                                  "Heliostat Field",                          "?=1.99",                                                           "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "eta_map",                            "Field efficiency array",                                                                                                                  "",             "",                                  "Heliostat Field",                          "?",                                                                "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "flux_maps",                          "Flux map intensities",                                                                                                                    "",             "",                                  "Heliostat Field",                          "?",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "c_atm_0",                            "Attenuation coefficient 0",                                                                                                               "",             "",                                  "Heliostat Field",                          "?=0.006789",                                                       "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "c_atm_1",                            "Attenuation coefficient 1",                                                                                                               "",             "",                                  "Heliostat Field",                          "?=0.1046",                                                         "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "c_atm_2",                            "Attenuation coefficient 2",                                                                                                               "",             "",                                  "Heliostat Field",                          "?=-0.0107",                                                        "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "c_atm_3",                            "Attenuation coefficient 3",                                                                                                               "",             "",                                  "Heliostat Field",                          "?=0.002845",                                                       "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "water_usage_per_wash",               "Water usage per wash",                                                                                                                    "L/m2_aper",    "",                                  "Heliostat Field",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "washing_frequency",                  "Mirror washing frequency",                                                                                                                "none",         "",                                  "Heliostat Field",                          "*",                                                                "",              ""},


    //{ SSC_INPUT,     SSC_NUMBER, "sf_excess",                          "Heliostat field multiple",                                                                                                                "",             "",                                  "System Design",                            "?=1.0",                                                            "",              ""},
    //{ SSC_INPUT,     SSC_NUMBER, "calc_fluxmaps",                      "Include fluxmap calculations",                                                                                                            "",             "",                                  "Heliostat Field",                          "?=0",                                                              "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "tower_fixed_cost",                   "Tower fixed cost",                                                                                                                        "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "tower_exp",                          "Tower cost scaling exponent",                                                                                                             "",             "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "rec_ref_cost",                       "Receiver reference cost",                                                                                                                 "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "rec_ref_area",                       "Receiver reference area for cost scale",                                                                                                  "",             "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "rec_cost_exp",                       "Receiver cost scaling exponent",                                                                                                          "",             "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "site_spec_cost",                     "Site improvement cost",                                                                                                                   "$/m2",         "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "heliostat_spec_cost",                "Heliostat field cost",                                                                                                                    "$/m2",         "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "plant_spec_cost",                    "Power cycle specific cost",                                                                                                               "$/kWe",        "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "bop_spec_cost",                      "BOS specific cost",                                                                                                                       "$/kWe",        "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "tes_spec_cost",                      "Thermal energy storage cost",                                                                                                             "$/kWht",       "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "land_spec_cost",                     "Total land area cost",                                                                                                                    "$/acre",       "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "contingency_rate",                   "Contingency for cost overrun",                                                                                                            "%",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "sales_tax_rate",                     "Sales tax rate",                                                                                                                          "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "sales_tax_frac",                     "Percent of cost to which sales tax applies",                                                                                              "%",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "cost_sf_fixed",                      "Solar field fixed cost",                                                                                                                  "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "fossil_spec_cost",                   "Fossil system specific cost",                                                                                                             "$/kWe",        "",                                  "System Costs",                             "*",                                                                "",              ""},

    //other costs needed for optimization update
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.cost.epc.per_acre",           "EPC cost per acre",                                                                                                                       "$/acre",       "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.cost.epc.percent",            "EPC cost percent of direct",                                                                                                              "%",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.cost.epc.per_watt",           "EPC cost per watt",                                                                                                                       "$/W",          "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.cost.epc.fixed",              "EPC fixed",                                                                                                                               "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.cost.plm.percent",            "PLM cost percent of direct",                                                                                                              "%",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.cost.plm.per_watt",           "PLM cost per watt",                                                                                                                       "$/W",          "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.cost.plm.fixed",              "PLM fixed",                                                                                                                               "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.sf.fixed_land_area",          "Fixed land area",                                                                                                                         "acre",         "",                                  "Heliostat Field",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.sf.land_overhead_factor",     "Land overhead factor",                                                                                                                    "",             "",                                  "Heliostat Field",                          "*",                                                                "",              ""},

	// System Design
    { SSC_INPUT,     SSC_NUMBER, "T_htf_cold_des",                     "Cold HTF inlet temperature at design conditions",                                                                                         "C",            "",                                  "System Design",                            "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "T_htf_hot_des",                      "Hot HTF outlet temperature at design conditions",                                                                                         "C",            "",                                  "System Design",                            "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "P_ref",                              "Reference output electric power at design condition",                                                                                     "MW",           "",                                  "System Design",                            "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "design_eff",                         "Power cycle efficiency at design",                                                                                                        "none",         "",                                  "System Design",                            "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "tshours",                            "Equivalent full-load thermal storage hours",                                                                                              "hr",           "",                                  "System Design",                            "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "solarm",                             "Solar multiple",                                                                                                                          "-",            "",                                  "System Design",                            "*",                                                                "",              ""},

	// Receiver (type 222) parameters
    { SSC_INPUT,     SSC_NUMBER, "N_panels",                           "Number of individual panels on the receiver",                                                                                             "",             "",                                  "Tower and Receiver",                       "*",                                                                "INTEGER",       ""},
    { SSC_INPUT,     SSC_NUMBER, "d_tube_out",                         "The outer diameter of an individual receiver tube",                                                                                       "mm",           "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "th_tube",                            "The wall thickness of a single receiver tube",                                                                                            "mm",           "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "mat_tube",                           "Receiver tube material, 2=Stainless AISI316",                                                                                             "",             "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "rec_htf",                            "Receiver HTF, 17=Salt (60% NaNO3, 40% KNO3) 10=Salt (46.5% LiF 11.5% NaF 42% KF) 50=Lookup tables",                                       "",             "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "field_fl_props",                     "User defined field fluid property data",                                                                                                  "-",            "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "Flow_type",                          "Receiver flow pattern: see figure on SAM Receiver page",                                                                                  "",             "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "crossover_shift",                    "Number of panels shift in receiver crossover position",                                                                                   "",             "",                                  "Tower and Receiver",                       "?=0",                                                              "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "epsilon",                            "The emissivity of the receiver surface coating",                                                                                          "",             "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "hl_ffact",                           "The heat loss factor (thermal loss fudge factor)",                                                                                        "",             "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "f_rec_min",                          "Minimum receiver mass flow rate turn down fraction",                                                                                      "",             "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "rec_su_delay",                       "Fixed startup delay time for the receiver",                                                                                               "hr",           "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "rec_qf_delay",                       "Energy-based receiver startup delay (fraction of rated thermal power)",                                                                   "",             "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.rec.max_oper_frac",           "Maximum receiver mass flow rate fraction",                                                                                                "",             "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "eta_pump",                           "Receiver HTF pump efficiency",                                                                                                            "",             "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "piping_loss",                        "Thermal loss per meter of piping",                                                                                                        "Wt/m",         "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "piping_length_mult",                 "Piping length multiplier",                                                                                                                "",             "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "piping_length_const",                "Piping constant length",                                                                                                                  "m",            "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},

	// TES parameters - general
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.tes.init_hot_htf_percent",    "Initial fraction of available volume that is hot",                                                                                        "%",            "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "h_tank",                             "Total height of tank (height of HTF when tank is full)",                                                                                  "m",            "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "cold_tank_max_heat",                 "Rated heater capacity for cold tank heating",                                                                                             "MW",           "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "u_tank",                             "Loss coefficient from the tank",                                                                                                          "W/m2-K",       "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "tank_pairs",                         "Number of equivalent tank pairs",                                                                                                         "",             "",                                  "Thermal Storage",                          "*",                                                                "INTEGER",       ""},
    { SSC_INPUT,     SSC_NUMBER, "cold_tank_Thtr",                     "Minimum allowable cold tank HTF temperature",                                                                                             "C",            "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
	// TES parameters - 2 tank
    { SSC_INPUT,     SSC_NUMBER, "h_tank_min",                         "Minimum allowable HTF height in storage tank",                                                                                            "m",            "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "hot_tank_Thtr",                      "Minimum allowable hot tank HTF temperature",                                                                                              "C",            "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "hot_tank_max_heat",                  "Rated heater capacity for hot tank heating",                                                                                              "MW",           "",                                  "Thermal Storage",                          "*",                                                                "",              ""},

	// Power Cycle Inputs
    { SSC_INPUT,     SSC_NUMBER, "pc_config",                          "PC configuration 0=Steam Rankine (224), 1=user defined, 2=sCO2 Recompression (424)",                                                      "",             "",                                  "Power Cycle",                              "?=0",                                                              "INTEGER",       ""},
    { SSC_INPUT,     SSC_NUMBER, "pb_pump_coef",                       "Pumping power to move 1kg of HTF through PB loop",                                                                                        "kW/kg",        "",                                  "Power Cycle",                              "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "startup_time",                       "Time needed for power block startup",                                                                                                     "hr",           "",                                  "Power Cycle",                              "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "startup_frac",                       "Fraction of design thermal power needed for startup",                                                                                     "none",         "",                                  "Power Cycle",                              "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "cycle_max_frac",                     "Maximum turbine over design operation fraction",                                                                                          "",             "",                                  "Power Cycle",                              "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "cycle_cutoff_frac",                  "Minimum turbine operation fraction before shutdown",                                                                                      "",             "",                                  "Power Cycle",                              "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "q_sby_frac",                         "Fraction of thermal power required for standby",                                                                                          "",             "",                                  "Power Cycle",                              "*",                                                                "",              ""},

	// User Defined cycle
    { SSC_INPUT,     SSC_NUMBER, "ud_T_amb_des",                       "Ambient temperature at user-defined power cycle design point",                                                                            "C",            "",                                  "User Defined Power Cycle",                 "pc_config=1",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "ud_f_W_dot_cool_des",                "Percent of user-defined power cycle design gross output consumed by cooling",                                                             "%",            "",                                  "User Defined Power Cycle",                 "pc_config=1",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "ud_m_dot_water_cool_des",            "Mass flow rate of water required at user-defined power cycle design point",                                                               "kg/s",         "",                                  "User Defined Power Cycle",                 "pc_config=1",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "ud_T_htf_low",                       "Low level HTF inlet temperature for T_amb parametric",                                                                                    "C",            "",                                  "User Defined Power Cycle",                 "pc_config=1",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "ud_T_htf_high",                      "High level HTF inlet temperature for T_amb parametric",                                                                                   "C",            "",                                  "User Defined Power Cycle",                 "pc_config=1",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "ud_T_amb_low",                       "Low level ambient temperature for HTF mass flow rate parametric",                                                                         "C",            "",                                  "User Defined Power Cycle",                 "pc_config=1",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "ud_T_amb_high",                      "High level ambient temperature for HTF mass flow rate parametric",                                                                        "C",            "",                                  "User Defined Power Cycle",                 "pc_config=1",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "ud_m_dot_htf_low",                   "Low level normalized HTF mass flow rate for T_HTF parametric",                                                                            "",             "",                                  "User Defined Power Cycle",                 "pc_config=1",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "ud_m_dot_htf_high",                  "High level normalized HTF mass flow rate for T_HTF parametric",                                                                           "",             "",                                  "User Defined Power Cycle",                 "pc_config=1",                                                      "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "ud_T_htf_ind_od",                    "Off design table of user-defined power cycle performance formed from parametric on T_htf_hot [C]",                                        "",             "",                                  "User Defined Power Cycle",                 "pc_config=1",                                                      "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "ud_T_amb_ind_od",                    "Off design table of user-defined power cycle performance formed from parametric on T_amb [C]",                                            "",             "",                                  "User Defined Power Cycle",                 "pc_config=1",                                                      "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "ud_m_dot_htf_ind_od",                "Off design table of user-defined power cycle performance formed from parametric on m_dot_htf [ND]",                                       "",             "",                                  "User Defined Power Cycle",                 "pc_config=1",                                                      "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "ud_ind_od",                          "Off design user-defined power cycle performance as function of T_htf, m_dot_htf [ND], and T_amb",                                         "",             "",                                  "User Defined Power Cycle",                 "pc_config=1",                                                      "",              ""},

// sCO2 Powerblock (type 424) inputs
    { SSC_INPUT,     SSC_NUMBER, "sco2_cycle_config",                  "SCO2 cycle configuration, 1=recompression, 2=partial cooling",                                                                            "",             "",                                  "SCO2 Cycle",                               "pc_config=2",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "eta_c",                              "Isentropic efficiency of compressor(s)",                                                                                                  "",             "",                                  "SCO2 Cycle",                               "pc_config=2",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "eta_t",                              "Isentropic efficiency of turbine",                                                                                                        "",             "",                                  "SCO2 Cycle",                               "pc_config=2",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "recup_eff_max",                      "Maximum recuperator effectiveness",                                                                                                       "",             "",                                  "SCO2 Cycle",                               "pc_config=2",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "P_high_limit",                       "Upper pressure limit in cycle",                                                                                                           "MPa",          "",                                  "SCO2 Cycle",                               "pc_config=2",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "deltaT_PHX",                         "Design temperature difference in PHX",                                                                                                    "C",            "",                                  "SCO2 Cycle",                               "pc_config=2",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "fan_power_perc_net",                 "Percent of net cycle output used for fan power at design",                                                                                "%",            "",                                  "SCO2 Cycle",                               "pc_config=2",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "sco2_T_amb_des",                     "Ambient temperature at design point",                                                                                                     "C",            "",                                  "SCO2 Cycle",                               "pc_config=2",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "sco2_T_approach",                    "Temperature difference between main compressor CO2 inlet and ambient air",                                                                "C",            "",                                  "SCO2 Cycle",                               "pc_config=2",                                                      "",              ""},

		// sCO2 Powerblock pre-process
    { SSC_INPUT,     SSC_NUMBER, "is_sco2_preprocess",                 "Is sco2 off-design performance preprocessed? 1=yes",                                                                                      "",             "",                                  "SCO2 Cycle",                               "?=0",                                                              "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "sco2ud_T_htf_cold_calc",             "HTF cold temperature from sCO2 cycle des, may be different than T_htf_cold_des",                                                          "C",            "",                                  "SCO2 Cycle",                               "is_sco2_preprocess=1",                                             "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "sco2ud_T_htf_low",                   "Low level HTF inlet temperature for T_amb parametric",                                                                                    "C",            "",                                  "SCO2 Cycle",                               "is_sco2_preprocess=1",                                             "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "sco2ud_T_htf_high",                  "High level HTF inlet temperature for T_amb parametric",                                                                                   "C",            "",                                  "SCO2 Cycle",                               "is_sco2_preprocess=1",                                             "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "sco2ud_T_amb_low",                   "Low level ambient temperature for HTF mass flow rate parametric",                                                                         "C",            "",                                  "SCO2 Cycle",                               "is_sco2_preprocess=1",                                             "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "sco2ud_T_amb_high",                  "High level ambient temperature for HTF mass flow rate parametric",                                                                        "C",            "",                                  "SCO2 Cycle",                               "is_sco2_preprocess=1",                                             "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "sco2ud_m_dot_htf_low",               "Low level normalized HTF mass flow rate for T_HTF parametric",                                                                            "",             "",                                  "SCO2 Cycle",                               "is_sco2_preprocess=1",                                             "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "sco2ud_m_dot_htf_high",              "High level normalized HTF mass flow rate for T_HTF parametric",                                                                           "",             "",                                  "SCO2 Cycle",                               "is_sco2_preprocess=1",                                             "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "sco2ud_T_htf_ind_od",                "Off design table of user-defined power cycle performance formed from parametric on T_htf_hot [C]",                                        "",             "",                                  "SCO2 Cycle",                               "is_sco2_preprocess=1",                                             "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "sco2ud_T_amb_ind_od",                "Off design table of user-defined power cycle performance formed from parametric on T_amb [C]",                                            "",             "",                                  "SCO2 Cycle",                               "is_sco2_preprocess=1",                                             "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "sco2ud_m_dot_htf_ind_od",            "Off design table of user-defined power cycle performance formed from parametric on m_dot_htf [ND]",                                       "",             "",                                  "SCO2 Cycle",                               "is_sco2_preprocess=1",                                             "",              ""},

    { SSC_INPUT,     SSC_NUMBER, "_sco2_P_high_limit",                 "Preprocess input: upper pressure limit",                                                                                                  "MPa",          "",                                  "SCO2 Cycle",                               "is_sco2_preprocess=1",                                             "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "_sco2_P_ref",                        "Preprocess input: gross power output",                                                                                                    "MWe",          "sco2_pc_pre",                       "SCO2 Cycle",                               "",                                                                 "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "_sco2_T_amb_des",                    "Preprocess input: design ambient temperature",                                                                                            "C",            "",                                  "SCO2 Cycle",                               "is_sco2_preprocess=1",                                             "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "_sco2_T_approach",                   "Preprocess input: compressor approach temperature",                                                                                       "C",            "",                                  "SCO2 Cycle",                               "is_sco2_preprocess=1",                                             "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "_sco2_T_htf_hot_des",                "Preprocess input: HTF hot temperature",                                                                                                   "C",            "",                                  "SCO2 Cycle",                               "is_sco2_preprocess=1",                                             "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "_sco2_deltaT_PHX",                   "Preprocess input: PHX approach temperature",                                                                                              "C",            "",                                  "SCO2 Cycle",                               "is_sco2_preprocess=1",                                             "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "_sco2_design_eff",                   "Preprocess input: cycle thermal efficiency",                                                                                              "",             "",                                  "SCO2 Cycle",                               "is_sco2_preprocess=1",                                             "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "_sco2_eta_c",                        "Preprocess input: compressor isentropic efficiency",                                                                                      "",             "",                                  "SCO2 Cycle",                               "is_sco2_preprocess=1",                                             "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "_sco2_eta_t",                        "Preprocess input: turbine isentropic efficiency",                                                                                         "",             "",                                  "SCO2 Cycle",                               "is_sco2_preprocess=1",                                             "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "_sco2_recup_eff_max",                "Preprocess input: max recuperator effectiveness",                                                                                         "",             "",                                  "SCO2 Cycle",                               "is_sco2_preprocess=1",                                             "",              ""},

	// System Control
    { SSC_INPUT,     SSC_NUMBER, "time_start",                         "Simulation start time",                                                                                                                   "s",            "",                                  "System Control",                           "?=0",                                                              "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "time_stop",                          "Simulation stop time",                                                                                                                    "s",            "",                                  "System Control",                           "?=31536000",                                                       "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "time_steps_per_hour",                "Number of simulation time steps per hour",                                                                                                "",             "",                                  "System Control",                           "?=-1",                                                             "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "vacuum_arrays",                      "Allocate arrays for only the required number of steps",                                                                                   "",             "",                                  "System Control",                           "?=0",                                                              "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "pb_fixed_par",                       "Fixed parasitic load - runs at all times",                                                                                                "MWe/MWcap",    "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "aux_par",                            "Aux heater, boiler parasitic",                                                                                                            "MWe/MWcap",    "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "aux_par_f",                          "Aux heater, boiler parasitic - multiplying fraction",                                                                                     "",             "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "aux_par_0",                          "Aux heater, boiler parasitic - constant coefficient",                                                                                     "",             "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "aux_par_1",                          "Aux heater, boiler parasitic - linear coefficient",                                                                                       "",             "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "aux_par_2",                          "Aux heater, boiler parasitic - quadratic coefficient",                                                                                    "",             "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "bop_par",                            "Balance of plant parasitic power fraction",                                                                                               "MWe/MWcap",    "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "bop_par_f",                          "Balance of plant parasitic power fraction - mult frac",                                                                                   "",             "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "bop_par_0",                          "Balance of plant parasitic power fraction - const coeff",                                                                                 "",             "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "bop_par_1",                          "Balance of plant parasitic power fraction - linear coeff",                                                                                "",             "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "bop_par_2",                          "Balance of plant parasitic power fraction - quadratic coeff",                                                                             "",             "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_ARRAY,  "f_turb_tou_periods",                 "Dispatch logic for turbine load fraction",                                                                                                "",             "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "weekday_schedule",                   "12x24 CSP operation Time-of-Use Weekday schedule",                                                                                        "",             "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "weekend_schedule",                   "12x24 CSP operation Time-of-Use Weekend schedule",                                                                                        "",             "",                                  "System Control",                           "*",                                                                "",              ""},

	// Financial inputs
    { SSC_INPUT,     SSC_MATRIX, "dispatch_sched_weekday",             "PPA pricing weekday schedule, 12x24",                                                                                                     "",             "",                                  "Time of Delivery Factors",                 "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "dispatch_sched_weekend",             "PPA pricing weekend schedule, 12x24",                                                                                                     "",             "",                                  "Time of Delivery Factors",                 "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dispatch_factor1",                   "Dispatch payment factor 1",                                                                                                               "",             "",                                  "Time of Delivery Factors",                 "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dispatch_factor2",                   "Dispatch payment factor 2",                                                                                                               "",             "",                                  "Time of Delivery Factors",                 "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dispatch_factor3",                   "Dispatch payment factor 3",                                                                                                               "",             "",                                  "Time of Delivery Factors",                 "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dispatch_factor4",                   "Dispatch payment factor 4",                                                                                                               "",             "",                                  "Time of Delivery Factors",                 "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dispatch_factor5",                   "Dispatch payment factor 5",                                                                                                               "",             "",                                  "Time of Delivery Factors",                 "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dispatch_factor6",                   "Dispatch payment factor 6",                                                                                                               "",             "",                                  "Time of Delivery Factors",                 "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dispatch_factor7",                   "Dispatch payment factor 7",                                                                                                               "",             "",                                  "Time of Delivery Factors",                 "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dispatch_factor8",                   "Dispatch payment factor 8",                                                                                                               "",             "",                                  "Time of Delivery Factors",                 "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dispatch_factor9",                   "Dispatch payment factor 9",                                                                                                               "",             "",                                  "Time of Delivery Factors",                 "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "is_dispatch_series",                 "Use time-series dispatch factors",                                                                                                        "",             "",                                  "System Control",                           "?=0",                                                              "",              ""},
    { SSC_INPUT,     SSC_ARRAY,  "dispatch_series",                    "Time series dispatch factors",                                                                                                            "",             "",                                  "System Control",                           "",                                                                 "",              ""},

// Inputs required for user defined SF performance
    { SSC_INPUT,     SSC_NUMBER, "A_sf_in",                            "Solar field area",                                                                                                                        "m^2",          "",                                  "Heliostat Field",                          "",                                                                 "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "A_sf",                               "Solar field area",                                                                                                                        "m^2",          "",                                  "",                                         "*",                                                                "",              ""},


	// optimized outputs updated depending on run type 
    { SSC_INPUT,     SSC_NUMBER, "rec_height",                         "Receiver height",                                                                                                                         "m",            "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
	{ SSC_INPUT,     SSC_NUMBER, "rec_width",                          "The receiver panel width",																												  "m",            "",                                  "Tower and Receiver",                       "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "h_tower",                            "Tower height",                                                                                                                            "m",            "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "N_hel",                              "Number of heliostats",                                                                                                                    "",             "",                                  "Heliostat Field",                          "",                                                                 "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "land_area_base",                     "Base land area occupied by heliostats",                                                                                                   "acre",         "",                                  "Heliostat Field",                          "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.total_land_area",        "Total land area",                                                                                                                         "acre",         "",                                  "System Costs",                             "*",                                                                "",              ""},

    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.site_improvements",      "Site improvement cost",                                                                                                                   "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.heliostats",             "Heliostat cost",                                                                                                                          "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.tower",                  "Tower cost",                                                                                                                              "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.receiver",               "Receiver cost",                                                                                                                           "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.storage",                "TES cost",                                                                                                                                "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.power_block",            "Power cycle cost",                                                                                                                        "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.rad_field",              "Radiative field cost"                                                                                                                     "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.rad_fluid",              "Radiative fluid cost"                                                                                                                     "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.rad_storage",            "Cold storage cost"                                                                                                                        "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.bop",                    "BOP cost",                                                                                                                                "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.fossil",                 "Fossil backup cost",                                                                                                                      "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "ui_direct_subtotal",                 "Direct capital precontingency cost",                                                                                                      "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.contingency",            "Contingency cost",                                                                                                                        "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "total_direct_cost",                  "Total direct cost",                                                                                                                       "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.epc.total",              "EPC and owner cost",                                                                                                                      "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.plm.total",              "Total land cost",                                                                                                                         "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.sales_tax.total",        "Sales tax cost",                                                                                                                          "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "total_indirect_cost",                "Total indirect cost",                                                                                                                     "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "total_installed_cost",               "Total installed cost",                                                                                                                    "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.installed_per_capacity", "Estimated installed cost per cap",                                                                                                        "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},

		// Construction financing inputs/outputs (SSC variable table from cmod_cb_construction_financing)
    { SSC_INPUT,     SSC_NUMBER, "const_per_interest_rate1",           "Interest rate, loan 1",                                                                                                                   "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_interest_rate2",           "Interest rate, loan 2",                                                                                                                   "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_interest_rate3",           "Interest rate, loan 3",                                                                                                                   "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_interest_rate4",           "Interest rate, loan 4",                                                                                                                   "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_interest_rate5",           "Interest rate, loan 5",                                                                                                                   "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_months1",                  "Months prior to operation, loan 1",                                                                                                       "",             "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_months2",                  "Months prior to operation, loan 2",                                                                                                       "",             "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_months3",                  "Months prior to operation, loan 3",                                                                                                       "",             "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_months4",                  "Months prior to operation, loan 4",                                                                                                       "",             "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_months5",                  "Months prior to operation, loan 5",                                                                                                       "",             "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_percent1",                 "Percent of total installed cost, loan 1",                                                                                                 "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_percent2",                 "Percent of total installed cost, loan 2",                                                                                                 "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_percent3",                 "Percent of total installed cost, loan 3",                                                                                                 "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_percent4",                 "Percent of total installed cost, loan 4",                                                                                                 "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_percent5",                 "Percent of total installed cost, loan 5",                                                                                                 "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_upfront_rate1",            "Upfront fee on principal, loan 1",                                                                                                        "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_upfront_rate2",            "Upfront fee on principal, loan 2",                                                                                                        "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_upfront_rate3",            "Upfront fee on principal, loan 3",                                                                                                        "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_upfront_rate4",            "Upfront fee on principal, loan 4",                                                                                                        "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_upfront_rate5",            "Upfront fee on principal, loan 5",                                                                                                        "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_principal1",               "Principal, loan 1",                                                                                                                       "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_principal2",               "Principal, loan 2",                                                                                                                       "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_principal3",               "Principal, loan 3",                                                                                                                       "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_principal4",               "Principal, loan 4",                                                                                                                       "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_principal5",               "Principal, loan 5",                                                                                                                       "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_interest1",                "Interest cost, loan 1",                                                                                                                   "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_interest2",                "Interest cost, loan 2",                                                                                                                   "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_interest3",                "Interest cost, loan 3",                                                                                                                   "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_interest4",                "Interest cost, loan 4",                                                                                                                   "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_interest5",                "Interest cost, loan 5",                                                                                                                   "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_total1",                   "Total financing cost, loan 1",                                                                                                            "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_total2",                   "Total financing cost, loan 2",                                                                                                            "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_total3",                   "Total financing cost, loan 3",                                                                                                            "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_total4",                   "Total financing cost, loan 4",                                                                                                            "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_total5",                   "Total financing cost, loan 5",                                                                                                            "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_percent_total",            "Total percent of installed costs, all loans",                                                                                             "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_principal_total",          "Total principal, all loans",                                                                                                              "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_interest_total",           "Total interest costs, all loans",                                                                                                         "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "construction_financing_cost",        "Total construction financing cost",                                                                                                       "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},



	// ****************************************************************************************************************************************
	// Outputs here:
	// ****************************************************************************************************************************************
		// Simulation outputs
    { SSC_OUTPUT,    SSC_ARRAY,  "time_hr",                            "Time at end of timestep",                                                                                                                 "hr",           "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "solzen",                             "Resource solar zenith",                                                                                                                   "deg",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "solaz",                              "Resource solar azimuth",                                                                                                                  "deg",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "beam",                               "Resource beam normal irradiance",                                                                                                         "W/m2",         "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "tdry",                               "Resource dry Bulb temperature",                                                                                                           "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "twet",                               "Resource wet Bulb temperature",                                                                                                           "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "rh",                                 "Resource relative humidity",                                                                                                              "%",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "wspd",                               "Resource wind velocity",                                                                                                                  "m/s",          "",                                  "",                                         "*",                                                                "",              ""},

    { SSC_OUTPUT,    SSC_ARRAY,  "q_sf_inc",                           "Field incident thermal power",                                                                                                            "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "eta_field",                          "Field optical efficiency",                                                                                                                "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "defocus",                            "Field optical focus fraction",                                                                                                            "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "sf_adjust_out",                      "Field availability adjustment factor",                                                                                                    "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_rec_inc",                      "Receiver incident thermal power",                                                                                                         "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "eta_therm",                          "Receiver thermal efficiency",                                                                                                             "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "Q_thermal",                          "Receiver thermal power to HTF less piping loss",                                                                                          "MWt",          "",                                  "",                                         "*",                                                                "",              ""},

    { SSC_OUTPUT,    SSC_ARRAY,  "m_dot_rec",                          "Receiver mass flow rate",                                                                                                                 "kg/s",         "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_startup",                          "Receiver startup thermal energy consumed",                                                                                                "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_rec_in",                           "Receiver HTF inlet temperature",                                                                                                          "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_rec_out",                          "Receiver HTF outlet temperature",                                                                                                         "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_piping_losses",                    "Receiver header/tower piping losses",                                                                                                     "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_thermal_loss",                     "Receiver convection and emission losses",                                                                                                 "MWt",          "",                                  "",                                         "*",                                                                "",              ""},

		// Power cycle outputs
    { SSC_OUTPUT,    SSC_ARRAY,  "eta",                                "PC efficiency, gross",                                                                                                                    "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_pb",                               "PC input energy",                                                                                                                         "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "m_dot_pc",                           "PC HTF mass flow rate",                                                                                                                   "kg/s",         "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_pc_startup",                       "PC startup thermal energy",                                                                                                               "MWht",         "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_pc_startup",                   "PC startup thermal power",                                                                                                                "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "P_cycle",                            "PC electrical power output, gross",                                                                                                       "MWe",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_pc_in",                            "PC HTF inlet temperature",                                                                                                                "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_pc_out",                           "PC HTF outlet temperature",                                                                                                               "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "m_dot_water_pc",                     "PC water consumption, makeup + cooling",                                                                                                  "kg/s",         "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_cond_out",                         "PC condenser water outlet temperature",                                                                                                   "C",            "",                                  "PC",                                       "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_cold",                             "Cold storage cold temperature",                                                                                                           "C",            "",                                  "PC",                                       "?",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "m_cold",                             "Cold storage cold tank mass",                                                                                                             "kg",           "",                                  "PC",                                       "?",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "m_warm",                             "Cold storage warm tank mass",                                                                                                             "kg",           "",                                  "PC",                                       "?",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_warm",                             "Cold storage warm tank temperature",                                                                                                      "C",            "",                                  "PC",                                       "?",                                                                "",              ""},

		// Thermal energy storage outputs
    { SSC_OUTPUT,    SSC_ARRAY,  "tank_losses",                        "TES thermal losses",                                                                                                                      "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_heater",                           "TES freeze protection power",                                                                                                             "MWe",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_tes_hot",                          "TES hot temperature",                                                                                                                     "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_tes_cold",                         "TES cold temperature",                                                                                                                    "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dc_tes",                           "TES discharge thermal power",                                                                                                             "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_ch_tes",                           "TES charge thermal power",                                                                                                                "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "e_ch_tes",                           "TES charge state",                                                                                                                        "MWht",         "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "m_dot_tes_dc",                       "TES discharge mass flow rate",                                                                                                            "kg/s",         "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "m_dot_tes_ch",                       "TES charge mass flow rate",                                                                                                               "kg/s",         "",                                  "",                                         "*",                                                                "",              ""},

		// Parasitics outputs
    { SSC_OUTPUT,    SSC_ARRAY,  "pparasi",                            "Parasitic power heliostat drives",                                                                                                        "MWe",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "P_tower_pump",                       "Parasitic power receiver/tower HTF pump",                                                                                                 "MWe",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "htf_pump_power",                     "Parasitic power TES and cycle HTF pump",                                                                                                  "MWe",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "P_cooling_tower_tot",                "Parasitic power condenser operation",                                                                                                     "MWe",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "P_fixed",                            "Parasitic power fixed load",                                                                                                              "MWe",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "P_plant_balance_tot",                "Parasitic power generation-dependent load",                                                                                               "MWe",          "",                                  "",                                         "*",                                                                "",              ""},

		// System outputs
    { SSC_OUTPUT,    SSC_ARRAY,  "P_out_net",                          "Total electric power to grid",                                                                                                            "MWe",          "",                                  "",                                         "*",                                                                "",              ""},

		// Controller outputs
    { SSC_OUTPUT,    SSC_ARRAY,  "tou_value",                          "CSP operating time-of-use value",                                                                                                         "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "pricing_mult",                       "PPA price multiplier",                                                                                                                    "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "n_op_modes",                         "Operating modes in reporting timestep",                                                                                                   "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "op_mode_1",                          "1st operating mode",                                                                                                                      "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "op_mode_2",                          "2nd operating mode, if applicable",                                                                                                       "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "op_mode_3",                          "3rd operating mode, if applicable",                                                                                                       "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "m_dot_balance",                      "Relative mass flow balance error",                                                                                                        "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_balance",                          "Relative energy balance error",                                                                                                           "",             "",                                  "",                                         "*",                                                                "",              ""},

		// These outputs correspond to the first csp-solver timestep in the reporting timestep.
		//     Subsequent csp-solver timesteps within the same reporting timestep are not tracked
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_pc_sb",                        "Thermal power for PC standby",                                                                                                            "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_pc_min",                       "Thermal power for PC min operation",                                                                                                      "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_pc_max",                       "Max thermal power to PC",                                                                                                                 "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_pc_target",                    "Target thermal power to PC",                                                                                                              "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "is_rec_su_allowed",                  "Is receiver startup allowed",                                                                                                             "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "is_pc_su_allowed",                   "Is power cycle startup allowed",                                                                                                          "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "is_pc_sb_allowed",                   "Is power cycle standby allowed",                                                                                                          "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_est_cr_su",                    "Estimated receiver startup thermal power",                                                                                                "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_est_cr_on",                    "Estimated receiver thermal power TO HTF",                                                                                                 "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_est_tes_dc",                   "Estimated max TES discharge thermal power",                                                                                               "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_est_tes_ch",                   "Estimated max TES charge thermal power",                                                                                                  "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "operating_modes_a",                  "First 3 operating modes tried",                                                                                                           "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "operating_modes_b",                  "Next 3 operating modes tried",                                                                                                            "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "operating_modes_c",                  "Final 3 operating modes tried",                                                                                                           "",             "",                                  "",                                         "*",                                                                "",              ""},

    { SSC_OUTPUT,    SSC_ARRAY,  "gen",                                "Total electric power to grid with available derate",                                                                                      "kWe",          "",                                  "",                                         "*",                                                                "",              ""},

    { SSC_OUTPUT,    SSC_MATRIX, "ud_T_htf_ind_od_out",                "T_htf_hot cycle off design",                                                                                                              "",             "",                                  "",                                         "?=[[0,1,2,3,4,5,6,7,8,9,10,11,12][0,1,2,3,4,5,6,7,8,9,10,11,12]]", "",              "COL_LABEL=UDPC_T_HTF_HOT,ROW_LABEL=NO_ROW_LABEL"},
    { SSC_OUTPUT,    SSC_MATRIX, "ud_T_amb_ind_od_out",                "T_amb cycle off design",                                                                                                                  "",             "",                                  "",                                         "?=[[0,1,2,3,4,5,6,7,8,9,10,11,12][0,1,2,3,4,5,6,7,8,9,10,11,12]]", "",              "COL_LABEL=UDPC_T_AMB,ROW_LABEL=NO_ROW_LABEL"},
    { SSC_OUTPUT,    SSC_MATRIX, "ud_m_dot_htf_ind_od_out",            "M_dot_htf cycle off design",                                                                                                              "",             "",                                  "",                                         "?=[[0,1,2,3,4,5,6,7,8,9,10,11,12][0,1,2,3,4,5,6,7,8,9,10,11,12]]", "",              "COL_LABEL=UDPC_M_DOT_HTF,ROW_LABEL=NO_ROW_LABEL"},

	// Annual single-value outputs
    { SSC_OUTPUT,    SSC_NUMBER, "annual_energy",                      "Annual total electric power to grid",                                                                                                     "kWhe",         "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "annual_W_cycle_gross",               "Electrical source - power cycle gross output",                                                                                            "kWhe",         "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "annual_W_cooling_tower",             "Total of condenser operation parasitics",                                                                                                 "kWhe",         "",                                  "PC",                                       "*",                                                                "",              ""},

    { SSC_OUTPUT,    SSC_NUMBER, "conversion_factor",                  "Gross to net conversion factor",                                                                                                          "%",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "capacity_factor",                    "Capacity factor",                                                                                                                         "%",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "kwh_per_kw",                         "First year kWh/kW",                                                                                                                       "kWh/kW",       "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "annual_total_water_use",             "Total annual water usage, cycle + mirror washing",                                                                                        "m3",           "",                                  "",                                         "*",                                                                "",              ""},

    { SSC_OUTPUT,    SSC_NUMBER, "disp_objective_ann",                 "Annual sum of dispatch objective function value",                                                                                         "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "disp_iter_ann",                      "Annual sum of dispatch solver iterations",                                                                                                "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "disp_presolve_nconstr_ann",          "Annual sum of dispatch problem constraint count",                                                                                         "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "disp_presolve_nvar_ann",             "Annual sum of dispatch problem variable count",                                                                                           "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "disp_solve_time_ann",                "Annual sum of dispatch solver time",                                                                                                      "",             "",                                  "",                                         "*",                                                                "",              ""},


    var_info_invalid };

class cm_gen3gas : public compute_module
{
public:

	cm_gen3gas()
	{
		add_var_info(_cm_vtab_gen3gas);
		add_var_info(vtab_adjustment_factors);
        add_var_info(vtab_sf_adjustment_factors);
	} 

	bool relay_message(string &msg, double percent)
	{
		log(msg);
		return update(msg, (float)percent);
	}

	void exec() throw(general_error)
	{
		// Weather reader
		C_csp_weatherreader weather_reader;
		if (is_assigned("solar_resource_file")){
			weather_reader.m_weather_data_provider = make_shared<weatherfile>(as_string("solar_resource_file"));
			if (weather_reader.m_weather_data_provider->has_message()) log(weather_reader.m_weather_data_provider->message(), SSC_WARNING);
		}
		if (is_assigned("solar_resource_data")){
			weather_reader.m_weather_data_provider = make_shared<weatherdata>(lookup("solar_resource_data"));
			if (weather_reader.m_weather_data_provider->has_message()) log(weather_reader.m_weather_data_provider->message(), SSC_WARNING);
		}

		weather_reader.m_trackmode = 0;
		weather_reader.m_tilt = 0.0;
		weather_reader.m_azimuth = 0.0;
		// Initialize to get weather file info
		weather_reader.init();
		if (weather_reader.has_error()) throw exec_error("gen3gas", weather_reader.get_error());

		// Get info from the weather reader initialization
		double site_elevation = weather_reader.ms_solved_params.m_elev;		//[m]


		int tes_type = 1;

		assign("q_design", as_number("P_ref") / as_number("design_eff") * as_number("solarm"));

		// Set up "cmod_solarpilot.cpp" conversions as necessary
		assign("helio_optical_error", (ssc_number_t)(as_number("helio_optical_error_mrad")*1.E-3));

		// Set 'n_flux_x' and 'n_flux_y' here, for now
		assign("n_flux_y", 1);
		int n_rec_panels = as_integer("N_panels");
		assign("n_flux_x", (ssc_number_t)max(12, n_rec_panels));

		// Calculate system capacity instead of pass in
		double system_capacity = as_double("P_ref") * as_double("gross_net_conversion_factor") *1.E3;		//[kWe]

		// 'sf_model_type'
		// 0 = design field and tower/receiver geometry
		// 1 = design field
		// 2 = user field, calculate performance
		// 3 = user performance maps vs solar position
		int field_model_type = as_integer("field_model_type");
		//	... require #3 for gen3 gas
		if (field_model_type != 3)
			throw exec_error("Gen3 gas model", "Field performance must be specified by the user for this model -- automatic calculation not currenty enabled.");

		// Run solarpilot right away to update values as needed
		solarpilot_invoke spi(this);
		util::matrix_t<double> mt_eta_map;
		util::matrix_t<double> mt_flux_maps;

		assign("is_optimize", 0);
		bool is_optimize = false;
		assign("calc_fluxmaps", 0);
		// The following optional inputs must be set here:
		assign("A_sf", as_number("A_sf_in"));

		if( tes_type != 1 )
		{
			throw exec_error("Gen3 gas model", "Thermocline thermal energy storage is not yet supported by the new CSP Solver and Dispatch Optimization models.\n");
		}

        
        // Set steps per hour
		C_csp_solver_gas::S_sim_setup sim_setup;
		sim_setup.m_sim_time_start = as_double("time_start");		//[s] time at beginning of first time step
		sim_setup.m_sim_time_end = as_double("time_stop");          //[s] time at end of last time step
        
        int steps_per_hour = (int)as_double("time_steps_per_hour");		//[-]

        //if the number of steps per hour is not provided (=-1), then assign it based on the weather file step
        if( steps_per_hour < 0 )
        {
            double sph_d = 3600. / weather_reader.m_weather_data_provider->step_sec();
            steps_per_hour = (int)( sph_d + 1.e-5 );
            if( (double)steps_per_hour != sph_d )
                throw spexception("The time step duration must be evenly divisible within an hour.");
        }

        size_t n_steps_fixed = (size_t)steps_per_hour * 8760;	//[-]
        if( as_boolean("vacuum_arrays") )
        {
            n_steps_fixed = steps_per_hour * (size_t)( (sim_setup.m_sim_time_end - sim_setup.m_sim_time_start)/3600. );
        }
        //int n_steps_fixed = (int)( (sim_setup.m_sim_time_end - sim_setup.m_sim_time_start) * steps_per_hour / 3600. ) ; 
		sim_setup.m_report_step = 3600.0 / (double)steps_per_hour;	//[s]





		// ***********************************************
		// ***********************************************
		// Power cycle
		// ***********************************************
		// ***********************************************
		C_csp_power_cycle * p_csp_power_cycle;
		// Steam Rankine and User Defined power cycle classes
		C_pc_Rankine_indirect_224 user_pc;
		// sCO2 power cycle class
		C_pc_sco2 sco2_pc;

		// Logic to choose between steam and sco2 power cycle 
		int pb_tech_type = as_integer("pc_config");

		if (pb_tech_type == 0)
		{
			throw exec_error("Gen3 gas model", "Steam Rankine cycle not allowed for this model. Please use sCO2 model or user-defined model.");
		}
		else if (pb_tech_type == 1)
		{
			C_pc_Rankine_indirect_224::S_params *pc = &user_pc.ms_params;
			pc->m_P_ref = as_double("P_ref");
			pc->m_eta_ref = as_double("design_eff");
			pc->m_T_htf_hot_ref = as_double("T_htf_hot_des");
			pc->m_T_htf_cold_ref = as_double("T_htf_cold_des");
			pc->m_cycle_max_frac = as_double("cycle_max_frac");
			pc->m_cycle_cutoff_frac = as_double("cycle_cutoff_frac");
			pc->m_q_sby_frac = as_double("q_sby_frac");
			pc->m_startup_time = as_double("startup_time");
			pc->m_startup_frac = as_double("startup_frac");
			pc->m_htf_pump_coef = as_double("pb_pump_coef");
			pc->m_pc_fl = as_integer("rec_htf");							// power cycle HTF is same as receiver HTF
			pc->m_pc_fl_props = as_matrix("field_fl_props");

			pc->m_is_user_defined_pc = true;

			// User-Defined Cycle Parameters
			pc->m_T_amb_des = as_double("ud_T_amb_des");	//[C]
			pc->m_W_dot_cooling_des = as_double("ud_f_W_dot_cool_des") / 100.0*as_double("P_ref");	//[MWe]
			pc->m_m_dot_water_des = as_double("ud_m_dot_water_cool_des");		//[kg/s]

			// Also need lower and upper levels for the 3 independent variables...
			pc->m_T_htf_low = as_double("ud_T_htf_low");			//[C]
			pc->m_T_htf_high = as_double("ud_T_htf_high");			//[C]
			pc->m_T_amb_low = as_double("ud_T_amb_low");			//[C]
			pc->m_T_amb_high = as_double("ud_T_amb_high");			//[C]
			pc->m_m_dot_htf_low = as_double("ud_m_dot_htf_low");	//[-]
			pc->m_m_dot_htf_high = as_double("ud_m_dot_htf_high");	//[-]

			// User-Defined Cycle Off-Design Tables 
			pc->mc_T_htf_ind = as_matrix("ud_T_htf_ind_od");
			pc->mc_T_amb_ind = as_matrix("ud_T_amb_ind_od");
			pc->mc_m_dot_htf_ind = as_matrix("ud_m_dot_htf_ind_od");
            pc->mc_combined_ind = as_matrix("ud_ind_od");

			// Set pointer to parent class
			p_csp_power_cycle = &user_pc;
		}
		else if (pb_tech_type == 2)
		{
			int is_sco2_preprocess = as_integer("is_sco2_preprocess");

			if (is_sco2_preprocess == 1)
			{
                double comp_var_tol = 1.E-6;    //[-]

				double _sco2_P_high_limit = as_double("_sco2_P_high_limit");   //[MPa]
				double P_high_limit = as_double("P_high_limit");			   //[MPa]
				if (!std::isfinite(_sco2_P_high_limit) || are_values_sig_different(_sco2_P_high_limit, P_high_limit, comp_var_tol) )
				{
					throw exec_error("gen3gas", util::format("The upper pressure limit used to generate"
						" the preprocessed sCO2 cycle data, %lg [MPa], is not equal to the input upper pressure limit %lg [MPa]",
						_sco2_P_high_limit, P_high_limit));
				}

				double _sco2_P_ref = as_double("_sco2_P_ref");
				double P_ref_input = as_double("P_ref");
				if (!std::isfinite(_sco2_P_ref) || are_values_sig_different(_sco2_P_ref, P_ref_input, comp_var_tol) )
				{
					throw exec_error("gen3gas", util::format("The cycle gross power used to generate"
						" the preprocessed sCO2 cycle data, %lg [MWe], is not equal to the input cycle gross power %lg [MWe]",
						_sco2_P_ref, P_ref_input));
				}

				double _sco2_T_amb_des = as_double("_sco2_T_amb_des");
				double T_amb_des_input = as_double("sco2_T_amb_des");
				if (!std::isfinite(_sco2_T_amb_des) || are_values_sig_different(_sco2_T_amb_des, T_amb_des_input, comp_var_tol) )
				{
					throw exec_error("gen3gas", util::format("The design ambient temperature used to generate"
						" the preprocessed sCO2 cycle data, %lg [C], is not equal to the input design ambient temperature %lg [C]",
						_sco2_T_amb_des, T_amb_des_input));
				}

				double _sco2_T_approach = as_double("_sco2_T_approach");
				double T_approach_input = as_double("sco2_T_approach");
				if (!std::isfinite(_sco2_T_approach) || are_values_sig_different(_sco2_T_approach, T_approach_input, comp_var_tol) )
				{
					throw exec_error("gen3gas", util::format("The compressor approach temperature used to generate"
						" the preprocessed sCO2 cycle data, %lg [C], is not equal to the input compressor approach temperature %lg [C]",
						_sco2_T_approach, T_approach_input));
				}

				double _sco2_T_htf_hot_des = as_double("_sco2_T_htf_hot_des");
				double T_htf_hot_des_input = as_double("T_htf_hot_des");
				if (!std::isfinite(_sco2_T_htf_hot_des) || are_values_sig_different(_sco2_T_htf_hot_des, T_htf_hot_des_input, comp_var_tol) )
				{
					throw exec_error("gen3gas", util::format("The HTF hot temperature uesd to generate"
						" the preprocessed sCO2 cycle data, %lg [C], is not equal to the input HTF hot temperature",
						_sco2_T_htf_hot_des, T_htf_hot_des_input));
				}

				double _sco2_deltaT_PHX = as_double("_sco2_deltaT_PHX");
				double deltaT_PHX_input = as_double("deltaT_PHX");
				if (!std::isfinite(_sco2_deltaT_PHX) || are_values_sig_different(_sco2_deltaT_PHX, deltaT_PHX_input, comp_var_tol) )
				{
					throw exec_error("gen3gas", util::format("The PHX approach temperature used to generate"
						" the preprocessed sCO2 cycle data, %lg [C], is not equal to the input PHX approach temperature",
						_sco2_deltaT_PHX, deltaT_PHX_input));
				}

				double _sco2_design_eff = as_double("_sco2_design_eff");
				double design_eff_input = as_double("design_eff");
				if ( !std::isfinite(_sco2_design_eff) || are_values_sig_different(_sco2_design_eff, design_eff_input, comp_var_tol) )
				{
					throw exec_error("gen3gas", util::format("The thermal efficiency used to generate"
						" the preprocessed sCO2 cycle data, %lg, is not equal to the input thermal efficiency",
						_sco2_design_eff, design_eff_input));
				}

				double _sco2_eta_c = as_double("_sco2_eta_c");
				double eta_c_input = as_double("eta_c");
				if (!std::isfinite(_sco2_eta_c) || are_values_sig_different(_sco2_eta_c, eta_c_input, comp_var_tol) )
				{
					throw exec_error("gen3gas", util::format("The compressor isentropic efficiency used to generate"
						" the preprocessed sCO2 cycle data, %lg, is not equal to the input compressor isentropic efficiency",
						_sco2_eta_c, eta_c_input));
				}

				double _sco2_eta_t = as_double("_sco2_eta_t");
				double eta_t_input = as_double("eta_t");
				if (!std::isfinite(_sco2_eta_t) || are_values_sig_different(_sco2_eta_t, eta_t_input, comp_var_tol) )
				{
					throw exec_error("gen3gas", util::format("The turbine isentropic efficiency used to generate"
						" the preprocessed sCO2 cycle data, %lg, is not equal to the input turbine isentropic efficiency",
						_sco2_eta_t, eta_t_input));
				}

				double _sco2_recup_eff_max = as_double("_sco2_recup_eff_max");
				double recup_eff_max = as_double("recup_eff_max");
				if (!std::isfinite(_sco2_recup_eff_max) || are_values_sig_different(_sco2_recup_eff_max, recup_eff_max, comp_var_tol) )
				{
					throw exec_error("gen3gas", util::format("The max recuperator effectiveness used to generate"
						" the preprocessed sCO2 cycle data, %lg, is not equal to the input max recuperator effectiveness",
						_sco2_recup_eff_max, recup_eff_max));
				}

				// ****************************************************
				// Setup UDPC model
				// ****************************************************
				C_pc_Rankine_indirect_224::S_params *pc = &user_pc.ms_params;
				pc->m_P_ref = as_double("P_ref");					//[MWe]
				pc->m_eta_ref = as_double("design_eff");			//[-]
				pc->m_T_htf_hot_ref = as_double("T_htf_hot_des");	//[C]
				
				double T_htf_cold_calc = as_double("sco2ud_T_htf_cold_calc");		//[C]
				double T_rec_htf_cold = as_double("T_htf_cold_des");				//[C]
				if (T_htf_cold_calc != T_rec_htf_cold)
				{
					assign("T_htf_cold_des", T_htf_cold_calc);							//[C]
					log(util::format("\nThe user input receiver design HTF cold temperature, %lg [C], was reset"
						" to the calculated sCO2 cycle HTF cold return temperature, %lg [C]\n", T_rec_htf_cold, T_htf_cold_calc), SSC_WARNING);
				}
				pc->m_T_htf_cold_ref = T_htf_cold_calc;		//[C]
				
				pc->m_cycle_max_frac = as_double("cycle_max_frac");		//[-]
				
				double cycle_cutoff_frac_sys = as_double("cycle_cutoff_frac");
				double cycle_cutoff_frac_sco2 = as_double("sco2ud_m_dot_htf_low");
				if (cycle_cutoff_frac_sys < cycle_cutoff_frac_sco2)
				{
					log(util::format("The user input cutoff fraction, %lg, was reset to the minimum allowable cutoff fraction"
						" for this sCO2 cycle off-design model, %lg\n", cycle_cutoff_frac_sys, cycle_cutoff_frac_sco2), SSC_WARNING);
					cycle_cutoff_frac_sys = cycle_cutoff_frac_sco2;
					assign("cycle_cutoff_frac", cycle_cutoff_frac_sys);
				}
				pc->m_cycle_cutoff_frac = cycle_cutoff_frac_sys;		//[-]
				
				pc->m_q_sby_frac = as_double("q_sby_frac");				//[-]
				pc->m_startup_time = as_double("startup_time");			//[hr]
				pc->m_startup_frac = as_double("startup_frac");			//[-]
				pc->m_htf_pump_coef = as_double("pb_pump_coef");		//[kWe/kg/s]
				pc->m_pc_fl = as_integer("rec_htf");					//[-] power cycle HTF is same as receiver HTF
				pc->m_pc_fl_props = as_matrix("field_fl_props");		//[-]

				// User-Defined Cycle Parameters
				pc->m_is_user_defined_pc = true;

				// User-Defined Cycle Parameters
				pc->m_T_amb_des = as_double("sco2_T_amb_des");			//[C]
				pc->m_W_dot_cooling_des = as_double("fan_power_perc_net") / 100.0*as_double("P_ref");	//[MWe]
				pc->m_m_dot_water_des = 0.0;		//[kg/s]

				// Also need lower and upper levels for the 3 independent variables...
				pc->m_T_htf_low = as_double("sco2ud_T_htf_low");			//[C]
				pc->m_T_htf_high = as_double("sco2ud_T_htf_high");			//[C]
				pc->m_T_amb_low = as_double("sco2ud_T_amb_low");			//[C]
				pc->m_T_amb_high = as_double("sco2ud_T_amb_high");			//[C]
				pc->m_m_dot_htf_low = as_double("sco2ud_m_dot_htf_low");	//[-]
				pc->m_m_dot_htf_high = as_double("sco2ud_m_dot_htf_high");	//[-]

				// User-Defined Cycle Off-Design Tables 
				pc->mc_T_htf_ind = as_matrix("sco2ud_T_htf_ind_od");
				pc->mc_T_amb_ind = as_matrix("sco2ud_T_amb_ind_od");
				pc->mc_m_dot_htf_ind = as_matrix("sco2ud_m_dot_htf_ind_od");

				p_csp_power_cycle = &user_pc;
			}
			else
			{
				// ****************************************
				// C_sco2_phx_air_cooler::S_des_par  User Defined Parameters
				// ****************************************
				C_sco2_phx_air_cooler::S_des_par sco2_rc_csp_par;
				// System Design Parameters
				sco2_rc_csp_par.m_hot_fl_code = as_integer("rec_htf");					//[-]
				sco2_rc_csp_par.mc_hot_fl_props = as_matrix("field_fl_props");			//[-]
				sco2_rc_csp_par.m_T_htf_hot_in = as_double("T_htf_hot_des") + 273.15;	//[K] Design HTF hot temp to power cycle
				sco2_rc_csp_par.m_phx_dt_hot_approach = as_double("deltaT_PHX");		//[K/C]
				sco2_rc_csp_par.m_T_amb_des = as_double("sco2_T_amb_des") + 273.15;		//[K] Design ambient temp, convert from C
				sco2_rc_csp_par.m_dt_mc_approach = as_double("sco2_T_approach");		//[K/C]
				sco2_rc_csp_par.m_elevation = site_elevation;							//[m]
				sco2_rc_csp_par.m_W_dot_net = as_double("P_ref")*1.E3;					//[kWe]

				sco2_rc_csp_par.m_cycle_config = as_integer("sco2_cycle_config");		//[-] 1 = recompression, 2 = partial cooling

				// Hardcode for now that design method iterates on UA_recup_total to hit target etas
				sco2_rc_csp_par.m_design_method = 1;
				sco2_rc_csp_par.m_eta_thermal = as_double("design_eff");				//[-]
				
				// Hardcode that recompression cycle is ok
				sco2_rc_csp_par.m_is_recomp_ok = 1;

				sco2_rc_csp_par.m_P_high_limit = as_double("P_high_limit")*1.E3;	//[kPa]
				sco2_rc_csp_par.m_fixed_P_mc_out = false;
				// Hardcode don't fix pressure ratio
				sco2_rc_csp_par.m_fixed_PR_HP_to_LP = false;				

				// ****************************************
				// ****************************************
				// C_sco2_phx_air_cooler::S_des_par  Hardcoded Parameters (for now...)
				// ****************************************
				// Cycle design parameters
				std::vector<double> DP_LT(2);
				/*(cold, hot) positive values are absolute [kPa], negative values are relative (-)*/
				DP_LT[0] = 0;
				DP_LT[1] = 0;
				/*(cold, hot) positive values are absolute [kPa], negative values are relative (-)*/
				std::vector<double> DP_HT(2);
				DP_HT[0] = 0;
				DP_HT[1] = 0;
				/*(cold, hot) positive values are absolute [kPa], negative values are relative (-)*/
				std::vector<double> DP_PC(2);
				DP_PC[0] = 0;
				DP_PC[1] = 0;
				/*(cold, hot) positive values are absolute [kPa], negative values are relative (-)*/
				std::vector<double> DP_PHX(2);
				DP_PHX[0] = 0;
				DP_PHX[1] = 0;
				sco2_rc_csp_par.m_DP_LT = DP_LT;
				sco2_rc_csp_par.m_DP_HT = DP_HT;
				sco2_rc_csp_par.m_DP_PC = DP_PC;
				sco2_rc_csp_par.m_DP_PHX = DP_PHX;
				sco2_rc_csp_par.m_N_sub_hxrs = 10;

				sco2_rc_csp_par.m_N_turbine = 3000.0;

				sco2_rc_csp_par.m_tol = 1.E-3;
				sco2_rc_csp_par.m_opt_tol = 1.E-3;
				
				// Cycle Design Parameters
                    // LTR thermal design
                sco2_rc_csp_par.m_LTR_target_code = 0;      // Optimize recuperator UA
                sco2_rc_csp_par.m_LTR_min_dT = std::numeric_limits<double>::quiet_NaN();    //[K] shouldn't be using this with target code = 1
				sco2_rc_csp_par.m_LTR_eff_max = as_double("recup_eff_max");		//[-]
                    // HTR thermal design
                sco2_rc_csp_par.m_HTR_target_code = 0;      // Optimize recuperator UA
                sco2_rc_csp_par.m_HTR_min_dT = std::numeric_limits<double>::quiet_NaN();    //[K] shouldn't be using this with target code = 1
				sco2_rc_csp_par.m_HTR_eff_max = as_double("recup_eff_max");		//[-]
                    //
				sco2_rc_csp_par.m_eta_mc = as_double("eta_c");					//[-]
				sco2_rc_csp_par.m_eta_rc = as_double("eta_c");					//[-]
				if (sco2_rc_csp_par.m_cycle_config == 2)
					sco2_rc_csp_par.m_eta_pc = as_double("eta_c");		   //[-]
				else
					sco2_rc_csp_par.m_eta_pc = sco2_rc_csp_par.m_eta_mc;
				sco2_rc_csp_par.m_eta_t = as_double("eta_t");					//[-]

				// PHX design parameters
				sco2_rc_csp_par.m_des_objective_type = 1;		//[-] Optimize design to maximize efficiency
				sco2_rc_csp_par.m_phx_dt_cold_approach = sco2_rc_csp_par.m_phx_dt_hot_approach;	//[K/C]

				// Air cooler parameters
				sco2_rc_csp_par.m_frac_fan_power = as_double("fan_power_perc_net") / 100.0;	//[-]
				sco2_rc_csp_par.m_deltaP_cooler_frac = 0.002;		//[-]

				sco2_pc.ms_params.ms_mc_sco2_recomp_params = sco2_rc_csp_par;

				bool is_preprocess_udpc = true;		// "is_preprocess_udpc"

				if (is_preprocess_udpc)
				{
					// For try/catch below
					int out_type = -1;
					std::string out_msg = "";

					//log("Calculating sCO2 design point...", SSC_WARNING);
					update("Calculating sCO2 design point...", 0.0);

					// Construction class and design system
                    C_sco2_phx_air_cooler c_sco2_csp;

					// Pass through callback function and pointer
                    c_sco2_csp.mf_callback_update = ssc_cmod_update;
                    c_sco2_csp.mp_mf_update = (void*)(this);

					try
					{
                        c_sco2_csp.design(sco2_rc_csp_par);
					}
					catch (C_csp_exception &csp_exception)
					{
						// Report warning before exiting with error
						while (c_sco2_csp.mc_messages.get_message(&out_type, &out_msg))
						{
							log(out_msg + "\n");
							log("\n");
						}

						throw exec_error("sco2_csp_system", csp_exception.m_error_message);
					}

					// Get sCO2 design outputs
					double T_htf_cold_calc = c_sco2_csp.get_design_solved()->ms_phx_des_solved.m_T_h_out;		//[K]
					log("sCO2 design point calculations complete.", SSC_WARNING);
					double T_rec_htf_cold = as_double("T_htf_cold_des");			//[C]
					assign("T_htf_cold_des", T_htf_cold_calc - 273.15);				//[C]
					log(util::format("\nThe user input receiver design HTF cold temperature, %lg [C], was reset"
						" to the calculated sCO2 cycle HTF cold return temperature, %lg [C]\n", T_rec_htf_cold, T_htf_cold_calc - 273.15), SSC_WARNING);
					update("Preprocessing cycle off-design...", 0.0);

					// Get user-defined power cycle parameters
					// HTF temperature parametric
					double T_htf_hot_low = c_sco2_csp.get_design_par()->m_T_htf_hot_in - 273.15 - 20.0;	//[C]
					double T_htf_hot_high = c_sco2_csp.get_design_par()->m_T_htf_hot_in - 273.15 + 15.0;	//[C]
					int n_T_htf_hot_in = 5;				//[-]

					// Ambient temperature parametric
					double T_amb_low = 0.0;				//[C]
					double T_amb_high = std::max(sco2_rc_csp_par.m_T_amb_des - 273.15 + 5.0, 45.0);			//[C]
					int n_T_amb_in = 10;				//[-]

					// HTF mass flow rate parametric
					double cycle_f_min = as_double("cycle_cutoff_frac");		//[-]
					bool is_des_rc = c_sco2_csp.get_design_solved()->ms_rc_cycle_solved.m_is_rc;	//[-]
					double sco2_f_min = 0.5;
					std::string cycle_type = "recompression";
					if (!is_des_rc)
					{
						sco2_f_min = 0.7;
						cycle_type = "simple";
					}
					if (cycle_f_min < sco2_f_min)
					{
						log(util::format("The user input cutoff fraction, %lg, was reset to the minimum allowable cutoff fraction"
							" for this sCO2 %s cycle off-design model, %lg\n", cycle_f_min, cycle_type.c_str(), sco2_f_min), SSC_WARNING);
						update("Preprocessing cycle off-design...", 0.0);
						cycle_f_min = sco2_f_min;
						assign("cycle_cutoff_frac", cycle_f_min);
					}
					double m_dot_htf_ND_low = min(0.95, cycle_f_min);	//[-]
					// Design is always = 1.0, so high needs to be a value > 1.0
					double m_dot_htf_ND_high = 1.05;			// max(1.05, as_double("cycle_max_frac"));		// + 0.01;		//[-]
					int n_m_dot_htf_ND_in = 10;

					util::matrix_t<double> T_htf_parametrics, T_amb_parametrics, m_dot_htf_ND_parametrics;

					try
					{
                        c_sco2_csp.generate_ud_pc_tables(T_htf_hot_low, T_htf_hot_high, n_T_htf_hot_in,
							T_amb_low, T_amb_high, n_T_amb_in,
							m_dot_htf_ND_low, m_dot_htf_ND_high, n_m_dot_htf_ND_in,
							T_htf_parametrics, T_amb_parametrics, m_dot_htf_ND_parametrics);
					}
					catch (C_csp_exception &csp_exception)
					{
						// Report warning before exiting with error
						while (c_sco2_csp.mc_messages.get_message(&out_type, &out_msg))
						{
							log(out_msg);
						}

						throw exec_error("sco2_csp_system", csp_exception.m_error_message);
					}

					size_t ncols = T_htf_parametrics.ncols();

					util::matrix_t<double> &p_udpc_T_htf_hot = allocate_matrix("ud_T_htf_ind_od_out", n_T_htf_hot_in, ncols);
					for (int i = 0; i < n_T_htf_hot_in; i++)
					{
						for (size_t j = 0; j < ncols; j++)
						{
							p_udpc_T_htf_hot(i, j) = (double)T_htf_parametrics(i, j);
						}
					}

					util::matrix_t<double> &p_udpc_T_amb = allocate_matrix("ud_T_amb_ind_od_out", n_T_amb_in, ncols);
					for (int i = 0; i < n_T_amb_in; i++)
					{
						for (size_t j = 0; j < ncols; j++)
						{
							p_udpc_T_amb(i, j) = (double)T_amb_parametrics(i, j);
						}
					}

					util::matrix_t<double> &p_udpc_m_dot_htf = allocate_matrix("ud_m_dot_htf_ind_od_out", n_m_dot_htf_ND_in, ncols);
					for (int i = 0; i < n_m_dot_htf_ND_in; i++)
					{
						for (size_t j = 0; j < ncols; j++)
						{
							p_udpc_m_dot_htf(i, j) = (double)m_dot_htf_ND_parametrics(i, j);
						}
					}

					log("sCO2 off-design performance calculations for lookup tables complete.", SSC_WARNING);
					update("sCO2 preprocess complete", 100.0);

					// ****************************************************
					// ****************************************************
					// Now, setup UDPC model
					// ****************************************************
					C_pc_Rankine_indirect_224::S_params *pc = &user_pc.ms_params;
					pc->m_P_ref = as_double("P_ref");
					pc->m_eta_ref = as_double("design_eff");
					pc->m_T_htf_hot_ref = as_double("T_htf_hot_des");
					pc->m_T_htf_cold_ref = c_sco2_csp.get_design_solved()->ms_phx_des_solved.m_T_h_out - 273.15;
					pc->m_cycle_max_frac = as_double("cycle_max_frac");
					pc->m_cycle_cutoff_frac = as_double("cycle_cutoff_frac");
					pc->m_q_sby_frac = as_double("q_sby_frac");
					pc->m_startup_time = as_double("startup_time");
					pc->m_startup_frac = as_double("startup_frac");
					pc->m_htf_pump_coef = as_double("pb_pump_coef");
					pc->m_pc_fl = as_integer("rec_htf");							// power cycle HTF is same as receiver HTF
					pc->m_pc_fl_props = as_matrix("field_fl_props");

					// User-Defined Cycle Parameters
					pc->m_is_user_defined_pc = true;

					pc->m_T_amb_des = c_sco2_csp.get_design_par()->m_T_amb_des - 273.15;	//[C]
					pc->m_W_dot_cooling_des = as_double("fan_power_perc_net") / 100.0*as_double("P_ref");	//[MWe]
					pc->m_m_dot_water_des = 0.0;		//[kg/s]

					// Also need lower and upper levels for the 3 independent variables...
					pc->m_T_htf_low = T_htf_hot_low;			//[C]
					pc->m_T_htf_high = T_htf_hot_high;			//[C]
					pc->m_T_amb_low = T_amb_low;				//[C]
					pc->m_T_amb_high = T_amb_high;				//[C]
					pc->m_m_dot_htf_low = m_dot_htf_ND_low;		//[-]
					pc->m_m_dot_htf_high = m_dot_htf_ND_high;	//[-]

					// User-Defined Cycle Off-Design Tables 
					pc->mc_T_htf_ind = T_htf_parametrics;
					pc->mc_T_amb_ind = T_amb_parametrics;
					pc->mc_m_dot_htf_ind = m_dot_htf_ND_parametrics;

					p_csp_power_cycle = &user_pc;
				}
				else
				{
					// ****************************************
					// ****************************************
					// C_sco2_phx_air_cooler::S_des_par   User Defined Parameters
					// ****************************************
					sco2_pc.ms_params.m_cycle_max_frac = as_double("cycle_max_frac");			//[-]
					sco2_pc.ms_params.m_cycle_cutoff_frac = as_double("cycle_cutoff_frac");		//[-]
					sco2_pc.ms_params.m_q_sby_frac = as_double("q_sby_frac");					//[-]
					sco2_pc.ms_params.m_startup_time = as_double("startup_time");				//[hr]
					sco2_pc.ms_params.m_startup_frac = as_double("startup_frac");				//[-]
					sco2_pc.ms_params.m_htf_pump_coef = as_double("pb_pump_coef");				//[kW/kg/s]

					p_csp_power_cycle = &sco2_pc;
				}
			}
		}
		else
		{
			std::string err_msg = util::format("The specified power cycle configuration, %d, does not exist. See SSC Input Table for help.\n", pb_tech_type);
			log(err_msg, SSC_WARNING);
			return;
		}

		// Set power cycle outputs common to all power cycle technologies
		p_csp_power_cycle->assign(C_pc_Rankine_indirect_224::E_ETA_THERMAL, allocate("eta", n_steps_fixed), n_steps_fixed);
		p_csp_power_cycle->assign(C_pc_Rankine_indirect_224::E_Q_DOT_HTF, allocate("q_pb", n_steps_fixed), n_steps_fixed);
		p_csp_power_cycle->assign(C_pc_Rankine_indirect_224::E_M_DOT_HTF, allocate("m_dot_pc", n_steps_fixed), n_steps_fixed);
		p_csp_power_cycle->assign(C_pc_Rankine_indirect_224::E_Q_DOT_STARTUP, allocate("q_dot_pc_startup", n_steps_fixed), n_steps_fixed);
		p_csp_power_cycle->assign(C_pc_Rankine_indirect_224::E_W_DOT, allocate("P_cycle", n_steps_fixed), n_steps_fixed);
		p_csp_power_cycle->assign(C_pc_Rankine_indirect_224::E_T_HTF_IN, allocate("T_pc_in", n_steps_fixed), n_steps_fixed);
		p_csp_power_cycle->assign(C_pc_Rankine_indirect_224::E_T_HTF_OUT, allocate("T_pc_out", n_steps_fixed), n_steps_fixed);
		p_csp_power_cycle->assign(C_pc_Rankine_indirect_224::E_M_DOT_WATER, allocate("m_dot_water_pc", n_steps_fixed), n_steps_fixed);
		p_csp_power_cycle->assign(C_pc_Rankine_indirect_224::E_T_COND_OUT, allocate("T_cond_out", n_steps_fixed), n_steps_fixed);

		//heliostat field class
		C_pt_sf_perf_interp heliostatfield;

		heliostatfield.ms_params.m_p_start = as_double("p_start");		//[kWe-hr] Heliostat startup energy
		heliostatfield.ms_params.m_p_track = as_double("p_track");		//[kWe] Heliostat tracking power
		heliostatfield.ms_params.m_hel_stow_deploy = as_double("hel_stow_deploy");	// N/A
		heliostatfield.ms_params.m_v_wind_max = as_double("v_wind_max");			// N/A
		heliostatfield.ms_params.m_n_flux_x = (int) as_double("n_flux_x");		// sp match
		heliostatfield.ms_params.m_n_flux_y = (int) as_double("n_flux_y");		// sp match

		heliostatfield.ms_params.m_eta_map = as_matrix("eta_map");
		heliostatfield.ms_params.m_flux_maps = as_matrix("flux_maps");
		heliostatfield.ms_params.m_N_hel = as_integer("N_hel");
		heliostatfield.ms_params.m_A_sf = as_double("A_sf");		//[m2]


        //Load the solar field adjustment factors
        sf_adjustment_factors sf_haf(this);
		size_t n_steps_full = weather_reader.m_weather_data_provider->nrecords(); //steps_per_hour * 8760;
		if (!sf_haf.setup((int)n_steps_full))
			throw exec_error("gen3gas", "failed to setup sf adjustment factors: " + sf_haf.error());
        //allocate array to pass to tcs
        heliostatfield.ms_params.m_sf_adjust.resize( sf_haf.size() );
        for( int i=0; i<sf_haf.size(); i++)     
            heliostatfield.ms_params.m_sf_adjust.at(i) = sf_haf(i);

		// Set callback information
		heliostatfield.mf_callback = ssc_cmod_solarpilot_callback;
		heliostatfield.m_cdata = (void*)this;

		// Try running pt heliostat init() call just for funsies
			// What happens when no callback to reference?
		//heliostatfield.init();


		//// *********************************************************
		////      Now set gas receiver parameters
		//// *********************************************************
		double H_rec = as_double("rec_height");
		double W_rec = as_double("rec_width");
		double rec_aspect = as_double("rec_aspect");

		double A_rec = std::numeric_limits<double>::quiet_NaN();

		A_rec = H_rec * W_rec * as_integer("N_panels");

		C_gen3gas_receiver receiver;
		receiver.m_n_panels = as_integer("N_panels");
		receiver.m_d_rec = W_rec;
		receiver.m_h_rec = H_rec;
		receiver.m_h_tower = as_double("h_tower");
		receiver.m_od_tube = as_double("d_tube_out");
		receiver.m_th_tube = as_double("th_tube");
		receiver.m_mat_tube = as_integer("mat_tube");
		receiver.m_field_fl = as_integer("rec_htf");
		receiver.m_field_fl_props = as_matrix("field_fl_props");
		receiver.m_flow_type = as_integer("Flow_type");
        receiver.m_crossover_shift = as_integer("crossover_shift");
		receiver.m_epsilon = as_double("epsilon");
		receiver.m_hl_ffact = as_double("hl_ffact");
		receiver.m_T_htf_hot_des = as_double("T_htf_hot_des");				//[C]
		receiver.m_T_htf_cold_des = as_double("T_htf_cold_des");			//[C]
		receiver.m_f_rec_min = as_double("f_rec_min");
		receiver.m_q_rec_des = as_double("P_ref")/as_double("design_eff")*as_double("solarm");
		receiver.m_rec_su_delay = as_double("rec_su_delay");
		receiver.m_rec_qf_delay = as_double("rec_qf_delay");
		receiver.m_m_dot_htf_max_frac = as_double("csp.pt.rec.max_oper_frac");
		receiver.m_A_sf = as_double("A_sf");

		// 8.10.2015 twn: add tower piping thermal losses to receiver performance
		receiver.m_pipe_loss_per_m = as_double("piping_loss");						//[Wt/m]
		receiver.m_pipe_length_add = as_double("piping_length_const");	//[m]
		receiver.m_pipe_length_mult = as_double("piping_length_mult");		//[-]

		receiver.m_n_flux_x = as_integer("n_flux_x");
		receiver.m_n_flux_y = as_integer("n_flux_y");

		receiver.m_T_salt_hot_target = as_double("T_htf_hot_des");
		receiver.m_eta_pump = as_double("eta_pump");
		receiver.m_night_recirc = 0;					// 8.15.15 twn: this is hardcoded for now - need to check that it is functioning correctly and reporting correct parasitics
		receiver.m_hel_stow_deploy = as_double("hel_stow_deploy");

		// Set parameters that were set with TCS defaults

		// Test mspt_receiver initialization
		//receiver.init();

		// Now try to instantiate mspt_collector_receiver
		C_csp_gas_collector_receiver collector_receiver(heliostatfield, receiver);
		// Then try init() call here, which should call inits from both classes
		//collector_receiver.init();

		// *******************************************************
		// *******************************************************
		// Set receiver outputs
		//float *p_q_thermal_copy = allocate("Q_thermal_123", n_steps_fixed);
		collector_receiver.mc_reported_outputs.assign(C_csp_gas_collector_receiver::E_FIELD_Q_DOT_INC, allocate("q_sf_inc", n_steps_fixed), n_steps_fixed);
		collector_receiver.mc_reported_outputs.assign(C_csp_gas_collector_receiver::E_FIELD_ETA_OPT, allocate("eta_field", n_steps_fixed), n_steps_fixed);
		collector_receiver.mc_reported_outputs.assign(C_csp_gas_collector_receiver::E_FIELD_ADJUST, allocate("sf_adjust_out", n_steps_fixed), n_steps_fixed);

		collector_receiver.mc_reported_outputs.assign(C_csp_gas_collector_receiver::E_Q_DOT_INC, allocate("q_dot_rec_inc", n_steps_fixed), n_steps_fixed);
		collector_receiver.mc_reported_outputs.assign(C_csp_gas_collector_receiver::E_ETA_THERMAL, allocate("eta_therm", n_steps_fixed), n_steps_fixed);
		collector_receiver.mc_reported_outputs.assign(C_csp_gas_collector_receiver::E_Q_DOT_THERMAL, allocate("Q_thermal", n_steps_fixed), n_steps_fixed);
		collector_receiver.mc_reported_outputs.assign(C_csp_gas_collector_receiver::E_M_DOT_HTF, allocate("m_dot_rec", n_steps_fixed), n_steps_fixed);
		collector_receiver.mc_reported_outputs.assign(C_csp_gas_collector_receiver::E_Q_DOT_STARTUP, allocate("q_startup", n_steps_fixed), n_steps_fixed);
		collector_receiver.mc_reported_outputs.assign(C_csp_gas_collector_receiver::E_T_HTF_IN, allocate("T_rec_in", n_steps_fixed), n_steps_fixed);
		collector_receiver.mc_reported_outputs.assign(C_csp_gas_collector_receiver::E_T_HTF_OUT, allocate("T_rec_out", n_steps_fixed), n_steps_fixed);
		collector_receiver.mc_reported_outputs.assign(C_csp_gas_collector_receiver::E_Q_DOT_PIPE_LOSS, allocate("q_piping_losses", n_steps_fixed), n_steps_fixed);
		collector_receiver.mc_reported_outputs.assign(C_csp_gas_collector_receiver::E_Q_DOT_LOSS, allocate("q_thermal_loss", n_steps_fixed), n_steps_fixed);


		// Thermal energy storage 
		C_csp_two_tank_tes storage;
		C_csp_two_tank_tes::S_params *tes = &storage.ms_params;
		tes->m_field_fl = as_integer("rec_htf");
		tes->m_field_fl_props = as_matrix("field_fl_props");
		tes->m_tes_fl = as_integer("rec_htf");
		tes->m_tes_fl_props = as_matrix("field_fl_props");
		tes->m_is_hx = false;									// MSPT assumes direct storage, so no user input required here: hardcode = false
		tes->m_W_dot_pc_design = as_double("P_ref");		//[MWe]
		tes->m_eta_pc = as_double("design_eff");				//[-]
		tes->m_solarm = as_double("solarm");
		tes->m_ts_hours = as_double("tshours");
		tes->m_h_tank = as_double("h_tank");
		tes->m_u_tank = as_double("u_tank");
		tes->m_tank_pairs = as_integer("tank_pairs");
		tes->m_hot_tank_Thtr = as_double("hot_tank_Thtr");
		tes->m_hot_tank_max_heat = as_double("hot_tank_max_heat");
		tes->m_cold_tank_Thtr = as_double("cold_tank_Thtr");
		tes->m_cold_tank_max_heat = as_double("cold_tank_max_heat");
		tes->m_dt_hot = 0.0;								// MSPT assumes direct storage, so no user input here: hardcode = 0.0
		tes->m_T_field_in_des = as_double("T_htf_cold_des");
		tes->m_T_field_out_des = as_double("T_htf_hot_des");
		tes->m_T_tank_hot_ini = as_double("T_htf_hot_des");
		tes->m_T_tank_cold_ini = as_double("T_htf_cold_des");
		tes->m_h_tank_min = as_double("h_tank_min");
		tes->m_f_V_hot_ini = as_double("csp.pt.tes.init_hot_htf_percent");
		tes->m_htf_pump_coef = as_double("pb_pump_coef");


		// TOU parameters
		C_csp_tou_block_schedules tou;
		C_csp_tou_block_schedules::S_params *tou_params = &tou.ms_params;
		tou_params->mc_csp_ops.mc_weekdays = as_matrix("weekday_schedule");
		tou_params->mc_csp_ops.mc_weekends = as_matrix("weekend_schedule");
		tou_params->mc_pricing.mc_weekdays = as_matrix("dispatch_sched_weekday");
		tou_params->mc_pricing.mc_weekends = as_matrix("dispatch_sched_weekend");
        /*tou.mc_dispatch_params.m_dispatch_optimize = as_boolean("is_dispatch");
        tou.mc_dispatch_params.m_is_write_ampl_dat = as_boolean("is_write_ampl_dat");
        tou.mc_dispatch_params.m_is_ampl_engine = as_boolean("is_ampl_engine");
        tou.mc_dispatch_params.m_ampl_data_dir = as_string("ampl_data_dir");
        tou.mc_dispatch_params.m_ampl_exec_call = as_string("ampl_exec_call");
		if( tou.mc_dispatch_params.m_dispatch_optimize )
		{
			tou.mc_dispatch_params.m_optimize_frequency = as_integer("disp_frequency");
            tou.mc_dispatch_params.m_disp_steps_per_hour = as_integer("disp_steps_per_hour");
			tou.mc_dispatch_params.m_optimize_horizon = as_integer("disp_horizon");
			tou.mc_dispatch_params.m_max_iterations = as_integer("disp_max_iter");
			tou.mc_dispatch_params.m_solver_timeout = as_double("disp_timeout");
			tou.mc_dispatch_params.m_mip_gap = as_double("disp_mip_gap");
			tou.mc_dispatch_params.m_presolve_type = as_integer("disp_spec_presolve");
			tou.mc_dispatch_params.m_bb_type = as_integer("disp_spec_bb");
			tou.mc_dispatch_params.m_disp_reporting = as_integer("disp_reporting");
			tou.mc_dispatch_params.m_scaling_type = as_integer("disp_spec_scaling");
			tou.mc_dispatch_params.m_disp_time_weighting = as_double("disp_time_weighting");
            tou.mc_dispatch_params.m_rsu_cost = as_double("disp_rsu_cost");
            tou.mc_dispatch_params.m_csu_cost = as_double("disp_csu_cost");
            tou.mc_dispatch_params.m_pen_delta_w = as_double("disp_pen_delta_w");
            tou.mc_dispatch_params.m_q_rec_standby = as_double("q_rec_standby");
			tou.mc_dispatch_params.m_w_rec_ht = as_double("q_rec_heattrace");

			if (as_boolean("is_wlim_series"))
			{
				size_t n_wlim_series = 0;
				ssc_number_t* wlim_series = as_array("wlim_series", &n_wlim_series);
				if (n_wlim_series != n_steps_full)
					throw exec_error("gen3gas", "Invalid net electricity generation limit series dimension. Matrix must have "+util::to_string((int)n_steps_full)+" rows.");
				for (size_t i = 0; i < n_steps_full; i++)
					tou.mc_dispatch_params.m_w_lim_full.at(i) = (double)wlim_series[i];
			}

	
		}*/
		tou.mc_dispatch_params.m_is_block_dispatch = ! tou.mc_dispatch_params.m_dispatch_optimize;      //mw
		tou.mc_dispatch_params.m_use_rule_1 = true;
		tou.mc_dispatch_params.m_standby_off_buffer = 2.0;
		tou.mc_dispatch_params.m_use_rule_2 = false;
		tou.mc_dispatch_params.m_q_dot_rec_des_mult = -1.23;
		tou.mc_dispatch_params.m_f_q_dot_pc_overwrite = -1.23;

        size_t n_f_turbine = 0;
		ssc_number_t *p_f_turbine = as_array("f_turb_tou_periods", &n_f_turbine);
		tou_params->mc_csp_ops.mvv_tou_arrays[C_block_schedule_csp_ops::TURB_FRAC].resize(n_f_turbine,0.0);
		//tou_params->mv_t_frac.resize(n_f_turbine, 0.0);
		for( size_t i = 0; i < n_f_turbine; i++ )
			tou_params->mc_csp_ops.mvv_tou_arrays[C_block_schedule_csp_ops::TURB_FRAC][i] = (double)p_f_turbine[i];

		bool is_timestep_input = (as_integer("ppa_multiplier_model") == 1);
		tou_params->mc_pricing.mv_is_diurnal = !(is_timestep_input);
		if (is_timestep_input)
		{
			size_t nmultipliers;
			ssc_number_t *multipliers = as_array("dispatch_factors_ts", &nmultipliers);
			tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE].resize(nmultipliers, 0.0);
			for (size_t ii = 0; ii < nmultipliers; ii++)
				tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][ii] = multipliers[ii];
		}
		else // standard diuranal input
		{
			tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE].resize(9, 0.0);
			tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][0] = as_double("dispatch_factor1");
			tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][1] = as_double("dispatch_factor2");
			tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][2] = as_double("dispatch_factor3");
			tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][3] = as_double("dispatch_factor4");
			tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][4] = as_double("dispatch_factor5");
			tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][5] = as_double("dispatch_factor6");
			tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][6] = as_double("dispatch_factor7");
			tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][7] = as_double("dispatch_factor8");
			tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][8] = as_double("dispatch_factor9");
		}

		// System parameters
		C_csp_solver_gas::S_csp_system_params system;
		system.m_pb_fixed_par = as_double("pb_fixed_par");
		system.m_bop_par = as_double("bop_par");
		system.m_bop_par_f = as_double("bop_par_f");
		system.m_bop_par_0 = as_double("bop_par_0");
		system.m_bop_par_1 = as_double("bop_par_1");
		system.m_bop_par_2 = as_double("bop_par_2");

  		// Instantiate Solver		
		C_csp_solver_gas csp_solver(weather_reader,
						collector_receiver, 
						*p_csp_power_cycle, 
						storage, 
						tou, 
						system,
						ssc_cmod_update,
						(void*)(this));


		// Set solver reporting outputs
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::TIME_FINAL, allocate("time_hr", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::ERR_M_DOT, allocate("m_dot_balance", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::ERR_Q_DOT, allocate("q_balance", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::N_OP_MODES, allocate("n_op_modes", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::OP_MODE_1, allocate("op_mode_1", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::OP_MODE_2, allocate("op_mode_2", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::OP_MODE_3, allocate("op_mode_3", n_steps_fixed), n_steps_fixed);


		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::TOU_PERIOD, allocate("tou_value", n_steps_fixed), n_steps_fixed);            
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::PRICING_MULT, allocate("pricing_mult", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::PC_Q_DOT_SB, allocate("q_dot_pc_sb", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::PC_Q_DOT_MIN, allocate("q_dot_pc_min", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::PC_Q_DOT_TARGET, allocate("q_dot_pc_target", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::PC_Q_DOT_MAX, allocate("q_dot_pc_max", n_steps_fixed), n_steps_fixed);
		
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::CTRL_IS_REC_SU, allocate("is_rec_su_allowed", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::CTRL_IS_PC_SU, allocate("is_pc_su_allowed", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::CTRL_IS_PC_SB, allocate("is_pc_sb_allowed", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::EST_Q_DOT_CR_SU, allocate("q_dot_est_cr_su", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::EST_Q_DOT_CR_ON, allocate("q_dot_est_cr_on", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::EST_Q_DOT_DC, allocate("q_dot_est_tes_dc", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::EST_Q_DOT_CH, allocate("q_dot_est_tes_ch", n_steps_fixed), n_steps_fixed);
		
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::CTRL_OP_MODE_SEQ_A, allocate("operating_modes_a", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::CTRL_OP_MODE_SEQ_B, allocate("operating_modes_b", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::CTRL_OP_MODE_SEQ_C, allocate("operating_modes_c", n_steps_fixed), n_steps_fixed);
		
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::SOLZEN, allocate("solzen", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::SOLAZ, allocate("solaz", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::BEAM, allocate("beam", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::TDRY, allocate("tdry", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::TWET, allocate("twet", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::RH, allocate("RH", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::WSPD, allocate("wspd", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::CR_DEFOCUS, allocate("defocus", n_steps_fixed), n_steps_fixed);

		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::TES_Q_DOT_LOSS, allocate("tank_losses", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::TES_W_DOT_HEATER, allocate("q_heater", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::TES_T_HOT, allocate("T_tes_hot", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::TES_T_COLD, allocate("T_tes_cold", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::TES_Q_DOT_DC, allocate("q_dc_tes", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::TES_Q_DOT_CH, allocate("q_ch_tes", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::TES_E_CH_STATE, allocate("e_ch_tes", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::TES_M_DOT_DC, allocate("m_dot_tes_dc", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::TES_M_DOT_CH, allocate("m_dot_tes_ch", n_steps_fixed), n_steps_fixed);

		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::COL_W_DOT_TRACK, allocate("pparasi", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::CR_W_DOT_PUMP, allocate("P_tower_pump", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::SYS_W_DOT_PUMP, allocate("htf_pump_power", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::PC_W_DOT_COOLING, allocate("P_cooling_tower_tot", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::SYS_W_DOT_FIXED, allocate("P_fixed", n_steps_fixed), n_steps_fixed);
		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::SYS_W_DOT_BOP, allocate("P_plant_balance_tot", n_steps_fixed), n_steps_fixed);

		csp_solver.mc_reported_outputs.assign(C_csp_solver_gas::C_solver_outputs::W_DOT_NET, allocate("P_out_net", n_steps_fixed), n_steps_fixed);



		update("Initialize Gen3 Gas model...", 0.0);

		int out_type = -1;
		std::string out_msg = "";
		try
		{
			// Initialize Solver
			csp_solver.init();
		}
		catch( C_csp_exception &csp_exception )
		{
			// Report warning before exiting with error
			while( csp_solver.mc_csp_messages.get_message(&out_type, &out_msg) )
			{
				log(out_msg, out_type);
			}

			throw exec_error("gen3gas", csp_exception.m_error_message);
		}

		// If no exception, then report messages
		while (csp_solver.mc_csp_messages.get_message(&out_type, &out_msg))
		{
			log(out_msg, out_type);
		}


        //if the pricing schedule is provided as hourly, overwrite the tou schedule
        if( as_boolean("is_dispatch_series") )
        {
            size_t n_dispatch_series;
            ssc_number_t *dispatch_series = as_array("dispatch_series", &n_dispatch_series);

            //resize the m_hr_tou array
            if( tou_params->mc_pricing.m_hr_tou != 0 )
                delete [] tou_params->mc_pricing.m_hr_tou;
            tou_params->mc_pricing.m_hr_tou = new double[n_steps_fixed];
            //set the tou period as unique for each time step
            for(size_t i=0; i<n_steps_fixed; i++)
                tou_params->mc_pricing.m_hr_tou[i] = i+1;
            //allocate reported arrays
            tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE].resize(n_steps_fixed);
            for( size_t i=0; i<n_steps_fixed; i++)
                tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][i] = dispatch_series[i];
        }

		update("Begin timeseries simulation...", 0.0);

		try
		{
			// Simulate !
			csp_solver.Ssimulate(sim_setup);
		}
		catch(C_csp_exception &csp_exception)
		{
			// Report warning before exiting with error
			while( csp_solver.mc_csp_messages.get_message(&out_type, &out_msg) )
			{
				log(out_msg);
			}

			throw exec_error("gen3gas", csp_exception.m_error_message);
		}

		// If no exception, then report messages
		while (csp_solver.mc_csp_messages.get_message(&out_type, &out_msg))
		{
			log(out_msg, out_type);
		}

		// ******* Re-calculate system costs here ************
		C_mspt_system_costs sys_costs;

		sys_costs.ms_par.A_sf_refl = as_double("A_sf");
		sys_costs.ms_par.site_improv_spec_cost = as_double("site_spec_cost");
		sys_costs.ms_par.heliostat_spec_cost = as_double("heliostat_spec_cost");
		sys_costs.ms_par.heliostat_fixed_cost = as_double("cost_sf_fixed");

		sys_costs.ms_par.h_tower = as_double("h_tower");
		sys_costs.ms_par.h_rec = H_rec;
		sys_costs.ms_par.h_helio = as_double("helio_height");
		sys_costs.ms_par.tower_fixed_cost = as_double("tower_fixed_cost");
		sys_costs.ms_par.tower_cost_scaling_exp = as_double("tower_exp");

		sys_costs.ms_par.A_rec = A_rec;
		sys_costs.ms_par.rec_ref_cost = as_double("rec_ref_cost");
		sys_costs.ms_par.A_rec_ref = as_double("rec_ref_area");
		sys_costs.ms_par.rec_cost_scaling_exp = as_double("rec_cost_exp");

		sys_costs.ms_par.Q_storage = as_double("P_ref") / as_double("design_eff")*as_double("tshours");
		sys_costs.ms_par.tes_spec_cost = as_double("tes_spec_cost");

		sys_costs.ms_par.W_dot_design = as_double("P_ref");
		sys_costs.ms_par.power_cycle_spec_cost = as_double("plant_spec_cost");

		sys_costs.ms_par.bop_spec_cost = as_double("bop_spec_cost");

		sys_costs.ms_par.fossil_backup_spec_cost = as_double("fossil_spec_cost");

		sys_costs.ms_par.contingency_rate = as_double("contingency_rate");

		//land area
		sys_costs.ms_par.total_land_area = as_double("land_area_base") * as_double("csp.pt.sf.land_overhead_factor") + as_double("csp.pt.sf.fixed_land_area")+ sys_costs.ms_par.radfield_area/4046.86 /*acres/m^2*/ ;
		assign("csp.pt.cost.total_land_area", (ssc_number_t)sys_costs.ms_par.total_land_area);

		sys_costs.ms_par.plant_net_capacity = system_capacity / 1000.0;			//[MWe], convert from kWe
		sys_costs.ms_par.EPC_land_spec_cost = as_double("csp.pt.cost.epc.per_acre");
		sys_costs.ms_par.EPC_land_perc_direct_cost = as_double("csp.pt.cost.epc.percent");
		sys_costs.ms_par.EPC_land_per_power_cost = as_double("csp.pt.cost.epc.per_watt");
		sys_costs.ms_par.EPC_land_fixed_cost = as_double("csp.pt.cost.epc.fixed");
		sys_costs.ms_par.total_land_spec_cost = as_double("land_spec_cost");
		sys_costs.ms_par.total_land_perc_direct_cost = as_double("csp.pt.cost.plm.percent");
		sys_costs.ms_par.total_land_per_power_cost = as_double("csp.pt.cost.plm.per_watt");
		sys_costs.ms_par.total_land_fixed_cost = as_double("csp.pt.cost.plm.fixed");
		sys_costs.ms_par.sales_tax_basis = as_double("sales_tax_frac");
		sys_costs.ms_par.sales_tax_rate = as_double("sales_tax_rate");

		try
		{
			sys_costs.calculate_costs();
		}
		catch (C_csp_exception &)
		{
			throw exec_error("MSPT system costs", util::format("System cost calculations failed. Check that all inputs are properly defined"));
		}

		// 1.5.2016 twn: financial model needs an updated total_installed_cost, remaining are for reporting only
		assign("total_installed_cost", (ssc_number_t)sys_costs.ms_out.total_installed_cost);

		assign("csp.pt.cost.site_improvements", (ssc_number_t)sys_costs.ms_out.site_improvement_cost);
		assign("csp.pt.cost.heliostats", (ssc_number_t)sys_costs.ms_out.heliostat_cost);
		assign("csp.pt.cost.tower", (ssc_number_t)sys_costs.ms_out.tower_cost);
		assign("csp.pt.cost.receiver", (ssc_number_t)sys_costs.ms_out.receiver_cost);
		assign("csp.pt.cost.storage", (ssc_number_t)sys_costs.ms_out.tes_cost);
		assign("csp.pt.cost.power_block", (ssc_number_t)sys_costs.ms_out.power_cycle_cost);
		
		assign("csp.pt.cost.bop", (ssc_number_t)sys_costs.ms_out.bop_cost);
		assign("csp.pt.cost.fossil", (ssc_number_t)sys_costs.ms_out.fossil_backup_cost);
		assign("ui_direct_subtotal", (ssc_number_t)sys_costs.ms_out.direct_capital_precontingency_cost);
		assign("csp.pt.cost.contingency", (ssc_number_t)sys_costs.ms_out.contingency_cost);
		assign("total_direct_cost", (ssc_number_t)sys_costs.ms_out.total_direct_cost);
		assign("csp.pt.cost.epc.total", (ssc_number_t)sys_costs.ms_out.epc_and_owner_cost);
		assign("csp.pt.cost.plm.total", (ssc_number_t)sys_costs.ms_out.total_land_cost);
		assign("csp.pt.cost.sales_tax.total", (ssc_number_t)sys_costs.ms_out.sales_tax_cost);
		assign("total_indirect_cost", (ssc_number_t)sys_costs.ms_out.total_indirect_cost);
		assign("csp.pt.cost.installed_per_capacity", (ssc_number_t)sys_costs.ms_out.estimated_installed_cost_per_cap);

		// Update construction financing costs, specifically, update: "construction_financing_cost"
		double const_per_interest_rate1 = as_double("const_per_interest_rate1");
		double const_per_interest_rate2 = as_double("const_per_interest_rate2");
		double const_per_interest_rate3 = as_double("const_per_interest_rate3");
		double const_per_interest_rate4 = as_double("const_per_interest_rate4");
		double const_per_interest_rate5 = as_double("const_per_interest_rate5");
		double const_per_months1 = as_double("const_per_months1");
		double const_per_months2 = as_double("const_per_months2");
		double const_per_months3 = as_double("const_per_months3");
		double const_per_months4 = as_double("const_per_months4");
		double const_per_months5 = as_double("const_per_months5");
		double const_per_percent1 = as_double("const_per_percent1");
		double const_per_percent2 = as_double("const_per_percent2");
		double const_per_percent3 = as_double("const_per_percent3");
		double const_per_percent4 = as_double("const_per_percent4");
		double const_per_percent5 = as_double("const_per_percent5");
		double const_per_upfront_rate1 = as_double("const_per_upfront_rate1");
		double const_per_upfront_rate2 = as_double("const_per_upfront_rate2");
		double const_per_upfront_rate3 = as_double("const_per_upfront_rate3");
		double const_per_upfront_rate4 = as_double("const_per_upfront_rate4");
		double const_per_upfront_rate5 = as_double("const_per_upfront_rate5");

		double const_per_principal1, const_per_principal2, const_per_principal3, const_per_principal4, const_per_principal5;
		double const_per_interest1, const_per_interest2, const_per_interest3, const_per_interest4, const_per_interest5;
		double const_per_total1, const_per_total2, const_per_total3, const_per_total4, const_per_total5;
		double const_per_percent_total, const_per_principal_total, const_per_interest_total, construction_financing_cost;

		const_per_principal1 = const_per_principal2 = const_per_principal3 = const_per_principal4 = const_per_principal5 =
			const_per_interest1 = const_per_interest2 = const_per_interest3 = const_per_interest4 = const_per_interest5 =
			const_per_total1 = const_per_total2 = const_per_total3 = const_per_total4 = const_per_total5 =
			const_per_percent_total = const_per_principal_total = const_per_interest_total = construction_financing_cost =
			std::numeric_limits<double>::quiet_NaN();

		N_financial_parameters::construction_financing_total_cost(sys_costs.ms_out.total_installed_cost,
			const_per_interest_rate1, const_per_interest_rate2, const_per_interest_rate3, const_per_interest_rate4, const_per_interest_rate5,
			const_per_months1, const_per_months2, const_per_months3, const_per_months4, const_per_months5,
			const_per_percent1, const_per_percent2, const_per_percent3, const_per_percent4, const_per_percent5,
			const_per_upfront_rate1, const_per_upfront_rate2, const_per_upfront_rate3, const_per_upfront_rate4, const_per_upfront_rate5,
			const_per_principal1, const_per_principal2, const_per_principal3, const_per_principal4, const_per_principal5,
			const_per_interest1, const_per_interest2, const_per_interest3, const_per_interest4, const_per_interest5,
			const_per_total1, const_per_total2, const_per_total3, const_per_total4, const_per_total5,
			const_per_percent_total, const_per_principal_total, const_per_interest_total, construction_financing_cost);

		assign("const_per_principal1", (ssc_number_t)const_per_principal1);
		assign("const_per_principal2", (ssc_number_t)const_per_principal2);
		assign("const_per_principal3", (ssc_number_t)const_per_principal3);
		assign("const_per_principal4", (ssc_number_t)const_per_principal4);
		assign("const_per_principal5", (ssc_number_t)const_per_principal5);
		assign("const_per_interest1", (ssc_number_t)const_per_interest1);
		assign("const_per_interest2", (ssc_number_t)const_per_interest2);
		assign("const_per_interest3", (ssc_number_t)const_per_interest3);
		assign("const_per_interest4", (ssc_number_t)const_per_interest4);
		assign("const_per_interest5", (ssc_number_t)const_per_interest5);
		assign("const_per_total1", (ssc_number_t)const_per_total1);
		assign("const_per_total2", (ssc_number_t)const_per_total2);
		assign("const_per_total3", (ssc_number_t)const_per_total3);
		assign("const_per_total4", (ssc_number_t)const_per_total4);
		assign("const_per_total5", (ssc_number_t)const_per_total5);
		assign("const_per_percent_total", (ssc_number_t)const_per_percent_total);
		assign("const_per_principal_total", (ssc_number_t)const_per_principal_total);
		assign("const_per_interest_total", (ssc_number_t)const_per_interest_total);
		assign("construction_financing_cost", (ssc_number_t)construction_financing_cost);

		// Do unit post-processing here
		double *p_q_pc_startup = allocate("q_pc_startup", n_steps_fixed);
		size_t count_pc_su = 0;
		ssc_number_t *p_q_dot_pc_startup = as_array("q_dot_pc_startup", &count_pc_su);
		if( count_pc_su != n_steps_fixed )
		{
			log("q_dot_pc_startup array is a different length than 'n_steps_fixed'.", SSC_WARNING);
			return;
		}
		for( size_t i = 0; i < n_steps_fixed; i++ )
		{
			p_q_pc_startup[i] = (float)(p_q_dot_pc_startup[i] * (sim_setup.m_report_step / 3600.0));	//[MWh]
		}

		// Convert mass flow rates from [kg/hr] to [kg/s]
		size_t count_m_dot_pc, count_m_dot_rec, count_m_dot_water_pc, count_m_dot_tes_dc, count_m_dot_tes_ch;
		count_m_dot_pc = count_m_dot_rec = count_m_dot_water_pc = count_m_dot_tes_dc = count_m_dot_tes_ch = 0;
		ssc_number_t *p_m_dot_rec = as_array("m_dot_rec", &count_m_dot_rec);
		ssc_number_t *p_m_dot_pc = as_array("m_dot_pc", &count_m_dot_pc);
		ssc_number_t *p_m_dot_water_pc = as_array("m_dot_water_pc", &count_m_dot_water_pc);
		ssc_number_t *p_m_dot_tes_dc = as_array("m_dot_tes_dc", &count_m_dot_tes_dc);
		ssc_number_t *p_m_dot_tes_ch = as_array("m_dot_tes_ch", &count_m_dot_tes_ch);
		if (count_m_dot_rec != n_steps_fixed || count_m_dot_pc != n_steps_fixed || count_m_dot_water_pc != n_steps_fixed
			|| count_m_dot_tes_dc != n_steps_fixed || count_m_dot_tes_ch != n_steps_fixed)
		{
			log("At least one m_dot array is a different length than 'n_steps_fixed'.", SSC_WARNING);
			return;
		}
		for (size_t i = 0; i < n_steps_fixed; i++)
		{
			p_m_dot_rec[i] = (ssc_number_t)(p_m_dot_rec[i] / 3600.0);	//[kg/s] convert from kg/hr
			p_m_dot_pc[i] = (ssc_number_t)(p_m_dot_pc[i] / 3600.0);		//[kg/s] convert from kg/hr
			p_m_dot_water_pc[i] = (ssc_number_t)(p_m_dot_water_pc[i] / 3600.0);	//[kg/s] convert from kg/hr
			p_m_dot_tes_dc[i] = (ssc_number_t)(p_m_dot_tes_dc[i] / 3600.0);		//[kg/s] convert from kg/hr
			p_m_dot_tes_ch[i] = (ssc_number_t)(p_m_dot_tes_ch[i] / 3600.0);		//[kg/s] convert from kg/hr
		}		

		size_t count;
		ssc_number_t *p_W_dot_net = as_array("P_out_net", &count);
		ssc_number_t *p_time_final_hr = as_array("time_hr", &count);

		// 'adjustment_factors' class stores factors in hourly array, so need to index as such
		adjustment_factors haf(this, "adjust");
		if( !haf.setup() )
			throw exec_error("gen3gas", "failed to setup adjustment factors: " + haf.error());


		ssc_number_t *p_gen = allocate("gen", count);
		for( size_t i = 0; i < count; i++ )
		{
			size_t hour = (size_t)ceil(p_time_final_hr[i]);
			p_gen[i] = (ssc_number_t)(p_W_dot_net[i] * 1.E3 * haf(hour));			//[kWe]
		}

		accumulate_annual_for_year("gen", "annual_energy", sim_setup.m_report_step / 3600.0, steps_per_hour, 1, n_steps_fixed/steps_per_hour);
		
		accumulate_annual_for_year("P_cycle", "annual_W_cycle_gross", 1000.0*sim_setup.m_report_step / 3600.0, steps_per_hour, 1, n_steps_fixed/steps_per_hour);		//[kWe-hr]
		accumulate_annual_for_year("P_cooling_tower_tot", "annual_W_cooling_tower", 1000.0*sim_setup.m_report_step / 3600.0, steps_per_hour, 1, n_steps_fixed / steps_per_hour);		//[kWe-hr]


		// Calculated Outputs
			// First, sum power cycle water consumption timeseries outputs
		accumulate_annual_for_year("m_dot_water_pc", "annual_total_water_use", sim_setup.m_report_step / 1000.0, steps_per_hour, 1, n_steps_fixed/steps_per_hour); //[m^3], convert from kg
			// Then, add water usage from mirror cleaning
		ssc_number_t V_water_cycle = as_number("annual_total_water_use");
		double V_water_mirrors = as_double("water_usage_per_wash") / 1000.0*as_double("A_sf")*as_double("washing_frequency");
		assign("annual_total_water_use", (ssc_number_t)(V_water_cycle + V_water_mirrors));

		ssc_number_t ae = as_number("annual_energy");
		ssc_number_t pg = as_number("annual_W_cycle_gross");
		ssc_number_t convfactor = (pg != 0) ? 100 * ae / pg : (ssc_number_t)0.0;
		assign("conversion_factor", convfactor);

		double kWh_per_kW = 0.0;
		double nameplate = system_capacity;		//[kWe]
		if(nameplate > 0.0)
			kWh_per_kW = ae / nameplate;

		assign("capacity_factor", (ssc_number_t)(kWh_per_kW / ((double)n_steps_fixed / (double)steps_per_hour)*100.));
		assign("kwh_per_kw", (ssc_number_t)kWh_per_kW);
		 
		//Single value outputs from radiative cooling system

	}
};

DEFINE_MODULE_ENTRY(gen3gas, "CSP Gen3 sCO2 gas power tower", 1)