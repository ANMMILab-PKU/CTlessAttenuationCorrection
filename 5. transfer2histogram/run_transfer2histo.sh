
./convert_Cdf_to_histo_ExternalMergeSort_OMP /share/home/xzhao/workDir/coin_sorter/sorter_20251021/histo/data_Cdf/0826/Coincidence_BrainPET_6p3p4_2DOI_Lu176_blankScan_20260408_0818_307_prompt_onlyTrue.Cdf /share/home/xzhao/workDir/coin_sorter/sorter_20251021/histo/data_histo/TBPET_1axial/0826/Coincidence_BrainPET_6p3p4_2DOI_Lu176_blankScan_20260408_0818_307_prompt_onlyTrue.histo

path1=/share/home/xzhao/workDir/coin_sorter/sorter_20251021/histo/data_histo/TBPET_1axial/0826/
file2=Coincidence_BrainPET_6p3p4_2DOI_Lu176_blankScan_20260408_0818_
file1=Coincidence_BrainPET_6p3p4_2DOI_Lu176_phantomScan_20260405_0818_
./histo_unify "${path1}${file1}202_prompt_onlyAll.histo" "${path1}${file1}202_prompt_onlyNotTrue.histo" "${path1}${file2}202_prompt_onlyTrue.histo"
./histo_unify "${path1}${file1}307_prompt_onlyAll.histo" "${path1}${file1}307_prompt_onlyNotTrue.histo" "${path1}${file2}307_prompt_onlyTrue.histo"



