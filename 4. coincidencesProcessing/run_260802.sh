
path2=/share/home/xzhao/workDir/coin_sorter/sorter_20251021/sorting/data_2608/
#file2=coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_IQP_w50mmBed_wSpheres_260729_0to900s_
file2=coincidence_TBPET_1axial_6p3p4_noDOI_Lu176_phantomScan_IQP_w50mmBed_wSpheres_260729_900to1800s_


 ./dat2Cdf_timeDiff_onlyTrue_sourceID0_disp "${path2}${file2}202_prompt.dat"
 ./dat2Cdf_timeDiff_onlyTrue_sourceID0_disp "${path2}${file2}307_prompt.dat"
 ./dat2Cdf_timeDiff_onlyTrue_sourceID1_disp "${path2}${file2}511_prompt.dat"

 mv "${path2}${file2}202_prompt.Cdf"   "${path2}${file2}202_prompt_onlyTrue.Cdf"
 mv "${path2}${file2}307_prompt.Cdf"   "${path2}${file2}307_prompt_onlyTrue.Cdf"
 mv "${path2}${file2}511_prompt.Cdf"   "${path2}${file2}511_prompt_onlyTrue.Cdf"
 mv "${path2}${file2}511_prompt.timeDiff"   "${path2}${file2}511_prompt_onlyTrue.timeDiff"


 ./dat2Cdf_timeDiff_onlyNotTrue_sourceID0_disp "${path2}${file2}202_prompt.dat"
 ./dat2Cdf_timeDiff_onlyNotTrue_sourceID0_disp "${path2}${file2}307_prompt.dat"
 ./dat2Cdf_timeDiff_onlyNotTrue_sourceID1_disp "${path2}${file2}511_prompt.dat"

 mv "${path2}${file2}202_prompt.Cdf"   "${path2}${file2}202_prompt_onlyNotTrue.Cdf"
 mv "${path2}${file2}307_prompt.Cdf"   "${path2}${file2}307_prompt_onlyNotTrue.Cdf"
 mv "${path2}${file2}511_prompt.Cdf"   "${path2}${file2}511_prompt_onlyNotTrue.Cdf"

 ./dat2Cdf_timeDiff_all_disp	   "${path2}${file2}202_prompt.dat"
 ./dat2Cdf_timeDiff_all_disp       "${path2}${file2}307_prompt.dat"
 ./dat2Cdf_timeDiff_all_disp       "${path2}${file2}511_prompt.dat"

 mv "${path2}${file2}202_prompt.Cdf"   "${path2}${file2}202_prompt_onlyAll.Cdf"
 mv "${path2}${file2}307_prompt.Cdf"   "${path2}${file2}307_prompt_onlyAll.Cdf"
 mv "${path2}${file2}511_prompt.Cdf"   "${path2}${file2}511_prompt_onlyAll.Cdf"
 mv "${path2}${file2}511_prompt.timeDiff"   "${path2}${file2}511_prompt_onlyAll.timeDiff"








