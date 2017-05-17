#!/bin/bash
RLE isData sample=merged_data15_periodD electron > RLE_merged_data15_elec_periodD.log 2>&1 &
RLE isData sample=merged_data15_periodE electron > RLE_merged_data15_elec_periodE.log 2>&1 &
RLE isData sample=merged_data15_periodF electron > RLE_merged_data15_elec_periodF.log 2>&1 &
RLE isData sample=merged_data15_periodG electron > RLE_merged_data15_elec_periodG.log 2>&1 &
RLE isData sample=merged_data15_periodH electron > RLE_merged_data15_elec_periodH.log 2>&1 &
RLE isData sample=merged_data15_periodJ electron > RLE_merged_data15_elec_periodJ.log 2>&1 &

RLE isData sample=merged_data16_periodA electron > RLE_merged_data16_periodA_elec.log 2>&1 &
RLE isData sample=merged_data16_periodB electron > RLE_merged_data16_periodB_elec.log 2>&1 &
RLE isData sample=merged_data16_periodC electron > RLE_merged_data16_periodC_elec.log 2>&1 &
RLE isData sample=merged_data16_periodD electron > RLE_merged_data16_periodD_elec.log 2>&1 &
RLE isData sample=merged_data16_periodE electron > RLE_merged_data16_periodE_elec.log 2>&1 &
RLE isData sample=merged_data16_periodF electron > RLE_merged_data16_periodF_elec.log 2>&1 &
RLE isData sample=merged_data16_periodG electron > RLE_merged_data16_periodG_elec.log 2>&1 &
RLE isData sample=merged_data16_periodI electron > RLE_merged_data16_periodI_elec.log 2>&1 &
RLE isData sample=merged_data16_periodK electron > RLE_merged_data16_periodK_elec.log 2>&1 &
RLE isData sample=merged_data16_periodL electron > RLE_merged_data16_periodL_elec.log 2>&1 &

RLE isData sample=merged_data15_periodD muon > RLE_merged_data15_muon_periodD.log 2>&1 &
RLE isData sample=merged_data15_periodE muon > RLE_merged_data15_muon_periodE.log 2>&1 &
RLE isData sample=merged_data15_periodF muon > RLE_merged_data15_muon_periodF.log 2>&1 &
RLE isData sample=merged_data15_periodG muon > RLE_merged_data15_muon_periodG.log 2>&1 &
RLE isData sample=merged_data15_periodH muon > RLE_merged_data15_muon_periodH.log 2>&1 &
RLE isData sample=merged_data15_periodJ muon > RLE_merged_data15_muon_periodJ.log 2>&1 &

RLE isData sample=merged_data16_periodA muon > RLE_merged_data16_periodA_muon.log 2>&1 &
RLE isData sample=merged_data16_periodB muon > RLE_merged_data16_periodB_muon.log 2>&1 &
RLE isData sample=merged_data16_periodC muon > RLE_merged_data16_periodC_muon.log 2>&1 &
RLE isData sample=merged_data16_periodD muon > RLE_merged_data16_periodD_muon.log 2>&1 &
RLE isData sample=merged_data16_periodE muon > RLE_merged_data16_periodE_muon.log 2>&1 &
RLE isData sample=merged_data16_periodF muon > RLE_merged_data16_periodF_muon.log 2>&1 &
RLE isData sample=merged_data16_periodG muon > RLE_merged_data16_periodG_muon.log 2>&1 &
RLE isData sample=merged_data16_periodI muon > RLE_merged_data16_periodI_muon.log 2>&1 &
RLE isData sample=merged_data16_periodK muon > RLE_merged_data16_periodK_muon.log 2>&1 &
RLE isData sample=merged_data16_periodL muon > RLE_merged_data16_periodL_muon.log 2>&1 &

RLE isMC Zee > RLE_Zee.log 2>&1 &
RLE isMC Zee_TandP_truth_match > RLE_Zee_TandP_truth_match.log 2>&1 &
RLE isMC Zee_truth_match > RLE_Zee_truth_match.log 2>&1 &
RLE isMC Zmumu > RLE_Zmumu.log 2>&1 &
RLE isMC Zmumu_TandP_truth_match > RLE_Zmumu_TandP_truth_match.log 2>&1 &
RLE isMC Zmumu_truth_match > RLE_Zmumu_truth_match.log 2>&1 &
RLE isMC ttbar electron > RLE_ttbar_elec.log 2>&1 &
RLE isMC ttbar muon > RLE_ttbar_muon.log 2>&1 &
# RLE isMC GG_ttn1 electron > RLE_GG_ttn1_elec.log 2>&1 &
# RLE isMC GG_ttn1 muon > RLE_GG_ttn1_muon.log 2>&1 &

# trigger uncertainties study. MC only
RLE isMC Zee tag_trigger_matched > RLE_Zee_tag_trigger_matched.log 2>&1 &
RLE isMC Zee trigger="dilepton_trigger" > RLE_Zee_dilepton_trigger.log 2>&1 &
RLE isMC Zee trigger="dilepton_trigger" tag_trigger_matched > RLE_Zee_dilepton_trigger_tag_trigger_matched.log 2>&1 &

RLE isMC Zmumu tag_trigger_matched > RLE_Zmumu_tag_trigger_matched.log 2>&1 &
RLE isMC Zmumu trigger="dilepton_trigger" > RLE_Zmumu_dilepton_trigger.log 2>&1 &
RLE isMC Zmumu trigger="dilepton_trigger" tag_trigger_matched > RLE_Zmumu_dilepton_trigger_tag_trigger_matched.log 2>&1 &

# For trigger uncertainties study. Data only.
# single lepton with tag trigger matched
RLE isData sample=merged_data15_periodD electron tag_trigger_matched > RLE_merged_data15_elec_periodD_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data15_periodE electron tag_trigger_matched > RLE_merged_data15_elec_periodE_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data15_periodF electron tag_trigger_matched > RLE_merged_data15_elec_periodF_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data15_periodG electron tag_trigger_matched > RLE_merged_data15_elec_periodG_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data15_periodH electron tag_trigger_matched > RLE_merged_data15_elec_periodH_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data15_periodJ electron tag_trigger_matched > RLE_merged_data15_elec_periodJ_tag_trigger_matched.log 2>&1 &

RLE isData sample=merged_data16_periodA electron tag_trigger_matched > RLE_merged_data16_periodA_elec_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodB electron tag_trigger_matched > RLE_merged_data16_periodB_elec_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodC electron tag_trigger_matched > RLE_merged_data16_periodC_elec_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodD electron tag_trigger_matched > RLE_merged_data16_periodD_elec_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodE electron tag_trigger_matched > RLE_merged_data16_periodE_elec_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodF electron tag_trigger_matched > RLE_merged_data16_periodF_elec_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodG electron tag_trigger_matched > RLE_merged_data16_periodG_elec_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodI electron tag_trigger_matched > RLE_merged_data16_periodI_elec_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodK electron tag_trigger_matched > RLE_merged_data16_periodK_elec_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodL electron tag_trigger_matched > RLE_merged_data16_periodL_elec_tag_trigger_matched.log 2>&1 &

RLE isData sample=merged_data15_periodD muon tag_trigger_matched > RLE_merged_data15_muon_periodD_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data15_periodE muon tag_trigger_matched > RLE_merged_data15_muon_periodE_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data15_periodF muon tag_trigger_matched > RLE_merged_data15_muon_periodF_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data15_periodG muon tag_trigger_matched > RLE_merged_data15_muon_periodG_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data15_periodH muon tag_trigger_matched > RLE_merged_data15_muon_periodH_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data15_periodJ muon tag_trigger_matched > RLE_merged_data15_muon_periodJ_tag_trigger_matched.log 2>&1 &

RLE isData sample=merged_data16_periodA muon tag_trigger_matched > RLE_merged_data16_periodA_muon_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodB muon tag_trigger_matched > RLE_merged_data16_periodB_muon_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodC muon tag_trigger_matched > RLE_merged_data16_periodC_muon_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodD muon tag_trigger_matched > RLE_merged_data16_periodD_muon_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodE muon tag_trigger_matched > RLE_merged_data16_periodE_muon_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodF muon tag_trigger_matched > RLE_merged_data16_periodF_muon_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodG muon tag_trigger_matched > RLE_merged_data16_periodG_muon_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodI muon tag_trigger_matched > RLE_merged_data16_periodI_muon_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodK muon tag_trigger_matched > RLE_merged_data16_periodK_muon_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodL muon tag_trigger_matched > RLE_merged_data16_periodL_muon_tag_trigger_matched.log 2>&1 &

# dilepton trigger
RLE isData sample=merged_data15_periodD electron trigger="dilepton_trigger" > RLE_merged_data15_elec_periodD_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data15_periodE electron trigger="dilepton_trigger" > RLE_merged_data15_elec_periodE_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data15_periodF electron trigger="dilepton_trigger" > RLE_merged_data15_elec_periodF_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data15_periodG electron trigger="dilepton_trigger" > RLE_merged_data15_elec_periodG_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data15_periodH electron trigger="dilepton_trigger" > RLE_merged_data15_elec_periodH_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data15_periodJ electron trigger="dilepton_trigger" > RLE_merged_data15_elec_periodJ_dilepton_trigger.log 2>&1 &

RLE isData sample=merged_data16_periodA electron trigger="dilepton_trigger" > RLE_merged_data16_periodA_elec_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data16_periodB electron trigger="dilepton_trigger" > RLE_merged_data16_periodB_elec_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data16_periodC electron trigger="dilepton_trigger" > RLE_merged_data16_periodC_elec_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data16_periodD electron trigger="dilepton_trigger" > RLE_merged_data16_periodD_elec_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data16_periodE electron trigger="dilepton_trigger" > RLE_merged_data16_periodE_elec_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data16_periodF electron trigger="dilepton_trigger" > RLE_merged_data16_periodF_elec_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data16_periodG electron trigger="dilepton_trigger" > RLE_merged_data16_periodG_elec_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data16_periodI electron trigger="dilepton_trigger" > RLE_merged_data16_periodI_elec_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data16_periodK electron trigger="dilepton_trigger" > RLE_merged_data16_periodK_elec_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data16_periodL electron trigger="dilepton_trigger" > RLE_merged_data16_periodL_elec_dilepton_trigger.log 2>&1 &

RLE isData sample=merged_data15_periodD muon trigger="dilepton_trigger" > RLE_merged_data15_muon_periodD_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data15_periodE muon trigger="dilepton_trigger" > RLE_merged_data15_muon_periodE_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data15_periodF muon trigger="dilepton_trigger" > RLE_merged_data15_muon_periodF_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data15_periodG muon trigger="dilepton_trigger" > RLE_merged_data15_muon_periodG_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data15_periodH muon trigger="dilepton_trigger" > RLE_merged_data15_muon_periodH_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data15_periodJ muon trigger="dilepton_trigger" > RLE_merged_data15_muon_periodJ_dilepton_trigger.log 2>&1 &

RLE isData sample=merged_data16_periodA muon trigger="dilepton_trigger" > RLE_merged_data16_periodA_muon_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data16_periodB muon trigger="dilepton_trigger" > RLE_merged_data16_periodB_muon_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data16_periodC muon trigger="dilepton_trigger" > RLE_merged_data16_periodC_muon_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data16_periodD muon trigger="dilepton_trigger" > RLE_merged_data16_periodD_muon_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data16_periodE muon trigger="dilepton_trigger" > RLE_merged_data16_periodE_muon_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data16_periodF muon trigger="dilepton_trigger" > RLE_merged_data16_periodF_muon_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data16_periodG muon trigger="dilepton_trigger" > RLE_merged_data16_periodG_muon_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data16_periodI muon trigger="dilepton_trigger" > RLE_merged_data16_periodI_muon_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data16_periodK muon trigger="dilepton_trigger" > RLE_merged_data16_periodK_muon_dilepton_trigger.log 2>&1 &
RLE isData sample=merged_data16_periodL muon trigger="dilepton_trigger" > RLE_merged_data16_periodL_muon_dilepton_trigger.log 2>&1 &

# dilepton trigger with tag trigger matched
RLE isData sample=merged_data15_periodD electron trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data15_elec_periodD_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data15_periodE electron trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data15_elec_periodE_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data15_periodF electron trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data15_elec_periodF_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data15_periodG electron trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data15_elec_periodG_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data15_periodH electron trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data15_elec_periodH_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data15_periodJ electron trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data15_elec_periodJ_dilepton_trigger_tag_trigger_matched.log 2>&1 &

RLE isData sample=merged_data16_periodA electron trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data16_periodA_elec_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodB electron trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data16_periodB_elec_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodC electron trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data16_periodC_elec_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodD electron trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data16_periodD_elec_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodE electron trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data16_periodE_elec_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodF electron trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data16_periodF_elec_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodG electron trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data16_periodG_elec_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodI electron trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data16_periodI_elec_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodK electron trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data16_periodK_elec_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodL electron trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data16_periodL_elec_dilepton_trigger_tag_trigger_matched.log 2>&1 &

RLE isData sample=merged_data15_periodD muon trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data15_muon_periodD_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data15_periodE muon trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data15_muon_periodE_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data15_periodF muon trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data15_muon_periodF_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data15_periodG muon trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data15_muon_periodG_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data15_periodH muon trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data15_muon_periodH_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data15_periodJ muon trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data15_muon_periodJ_dilepton_trigger_tag_trigger_matched.log 2>&1 &

RLE isData sample=merged_data16_periodA muon trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data16_periodA_muon_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodB muon trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data16_periodB_muon_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodC muon trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data16_periodC_muon_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodD muon trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data16_periodD_muon_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodE muon trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data16_periodE_muon_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodF muon trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data16_periodF_muon_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodG muon trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data16_periodG_muon_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodI muon trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data16_periodI_muon_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodK muon trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data16_periodK_muon_dilepton_trigger_tag_trigger_matched.log 2>&1 &
RLE isData sample=merged_data16_periodL muon trigger="dilepton_trigger" tag_trigger_matched > RLE_merged_data16_periodL_muon_dilepton_trigger_tag_trigger_matched.log 2>&1 &
