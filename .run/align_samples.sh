#for cline in 22Rv1 A549 HepG2 HeLa_Kyoto_MboI_G1_Snyc_STAG1_Depleted HeLa_Kyoto_MboI_G1_Snyc_STAG2_Depleted HeLa_MboI_G1_Snyc HL60S4_14_Pore HL60S4_5_Pore HL60S4 Nalm6 T47D;

# init file list
fns=""

###### ADDING SRA FILES ######
# For loop for samples with weird RE
#for cline in A549 HepG2;

# For loop for samples still downloading
#for cline in HeLa_Kyoto_MboI_G1_Snyc_STAG1_Depleted HeLa_Kyoto_MboI_G1_Snyc_STAG2_Depleted HeLa_MboI_G1_Snyc;

# For loop for samples with fully downloaded data
#for cline in 22Rv1 HL60S4 HL60S4_14_Pore HL60S4_5_Pore Nalm6;
#for cline in HL60S4_14_Pore Nalm6;
#for cline in HL60S4_14_Pore;
#for cline in 22Rv1;
#for cline in Nalm6;

#for cline in A549;
for cline in HepG2 A549 HL60S4 HL60S4_5_Pore HL60S4_14_Pore;
do
    # mapping
    fn="results/main/${cline}/hicpro_with_parallel/qsubs.started"

    # concat current file
    fns+="${fn} "
done

####### ADDING ENCODE FILES ######
#cline='T47D'
#libs=$(ls results/main/T47D/reads/T47D.ENC*.tsv | sed 's!.*/!!' | cut -f 2 -d ".")
#for lib in $libs;
#do
#
#    # mapping
#    #fn="results/main/${cline}/hicpro/bowtie_results/bwt2/${lib}/${lib}_1_hg38.bwt2merged.bam"


#    # running the whole pipeline
#    fn="results/main/${cline}/hicnv/${cline}.${lib}_hicnv/CNV_Estimation/${cline}.${lib}.cnv.bedGraph"
#    fns+=" ${fn}"
#    break
#done

echo "snakemake --profile workflow/profiles/pbs-torque/ $@ $fns"
snakemake --profile workflow/profiles/pbs-torque/ $@ $fns
