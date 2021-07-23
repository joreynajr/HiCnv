for r in R1 R2;
do 
     # smaller sam
     /mnt/BioApps/samtools-1.9/bin/samtools view -h \
            "results/main/LNCaP/4d_nucleome/4DNESNHN919R-B1-T1_${r}.bam" | \
            head -n 1000 > \
            "results/main/LNCaP/4d_nucleome/4DNESNHN919R-B1-T1_${r}.sorted.small.sam"

     # convert smaller sam to bam
     /mnt/BioApps/samtools-1.9/bin/samtools view -hb \
            "results/main/LNCaP/4d_nucleome/4DNESNHN919R-B1-T1_${r}.sorted.small.sam" > \
            "results/main/LNCaP/4d_nucleome/4DNESNHN919R-B1-T1_${r}.sorted.small.bam"

     # sort the smaller bam
     /mnt/BioApps/samtools-1.9/bin/samtools sort -n \
            "results/main/LNCaP/4d_nucleome/4DNESNHN919R-B1-T1_${r}.sorted.small.bam" > \
            "results/main/LNCaP/4d_nucleome/4DNESNHN919R-B1-T1_${r}.sorted.small.readSorted.bam"

     # sort the smaller bam
     /mnt/BioApps/samtools-1.9/bin/samtools sort -n \
            "results/main/LNCaP/4d_nucleome/4DNESNHN919R-B1-T1_${r}.sorted.small.bam" > \
            "results/main/LNCaP/4d_nucleome/4DNESNHN919R-B1-T1_${r}.sorted.small.readSorted.bam"

     # sort the smaller bam
     /mnt/BioApps/samtools-1.9/bin/samtools index \
            "results/main/LNCaP/4d_nucleome/4DNESNHN919R-B1-T1_${r}.sorted.small.readSorted.bam"

done

