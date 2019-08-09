#! /usr/bin/env bash
#File: assembly_pipeline.sh

##########
#Inputs
##########
config_file=assembly.config

##########
#functions
##########


function randomize_contigs {
	local read_array=($@)
        local read_header=SRR
	local path_data_fastq=/mnt/data/troyalty_data/sequencing/fastq_files/
        local read_tailer_forward=_dedupe_1.fastq
        local read_tailer_reverse=_dedupe_2.fastq
	local read_tailer_forward2=_1.fastq
	local read_tailer_reverse2=_2.fastq
	local read_tailer_fraction=_fraction
	local lines=
	local sample_reads=
	local total_reads=
	local seed=	

	for read in $read_array; do
		if [[ -f $path_data_fastq$read_header$read/$read_header$read$read_tailer_fraction$read_tailer_forward2 ]]; then
			rm $path_data_fastq$read_header$read/$read_header$read$read_tailer_fraction$read_tailer_forward2 $path_data_fastq$read_header$read/$read_header$read$read_tailer_fraction$read_tailer_reverse2
		fi

		seed=$(shuf -i 1-100000000 | head -n 1)
		fastq-sample $path_data_fastq$read_header$read/$read_header$read$read_tailer_forward $path_data_fastq$read_header$read/$read_header$read$read_tailer_reverse -p $fraction -o tmp -s $seed
		mv tmp.1.fastq $path_data_fastq$read_header$read/$read_header$read$read_tailer_fraction$read_tailer_forward2
		mv tmp.2.fastq $path_data_fastq$read_header$read/$read_header$read$read_tailer_fraction$read_tailer_reverse2
#		lines=($(wc -l $path_data_fastq$read_header$read/$read_header$read$read_tailer_forward))
#		total_reads=$(echo "$lines/4/1" | bc)
#		sample_reads=$(echo "$total_reads*$fraction/1" | bc)

#		seq 1 $total_reads | shuf | head -n $sample_reads > tmp
#		awk -F" " '{$1*=4;$1-=3;print}' tmp > tmp2
#		awk 'NR==FNR {a[$0]; a[$0+1]; a[$0+2]; a[$0+3]; next}; FNR in a' tmp2 $path_data_fastq$read_header$read/$read_header$read$read_tailer_forward > $path_data_fastq$read_header$read/$read_header$read$read_tailer_fraction$read_tailer_forward2
#		awk 'NR==FNR {a[$0]; a[$0+1]; a[$0+2]; a[$0+3]; next}; FNR in a' tmp2 $path_data_fastq$read_header$read/$read_header$read$read_tailer_reverse > $path_data_fastq$read_header$read/$read_header$read$read_tailer_fraction$read_tailer_reverse2
#		rm tmp*
		
	done
}



function remove_fastq_duplicates {
	local read_array=($@)
	local read_header=SRR
	local path_data_fastq=/mnt/data/troyalty_data/sequencing/fastq_files/
        local read_tailer_forward=_1.fastq
        local read_tailer_reverse=_2.fastq
	local dup=_dedupe

	for read in $read_array; do
		mv $path_data_fastq$read_header$read/$read_header$read$read_tailer_forward $path_data_fastq$read_header$read/tmp1
		mv $path_data_fastq$read_header$read/$read_header$read$read_tailer_reverse $path_data_fastq$read_header$read/tmp2
		perl ~/Software/prinseq-lite-0.20.4/prinseq-lite.pl -fastq $path_data_fastq$read_header$read/tmp1 -fastq2 $path_data_fastq$read_header$read/tmp2 -derep 1
		rm $path_data_fastq$read_header$read/*singletons*
		mv $path_data_fastq$read_header$read/*1_prinseq_good* $path_data_fastq$read_header$read/$read_header$read$dup$read_tailer_forward
		mv $path_data_fastq$read_header$read/*2_prinseq_good* $path_data_fastq$read_header$read/$read_header$read$dup$read_tailer_reverse
		rm $path_data_fastq$read_header$read/tmp*
		rm $path_data_fastq$read_header$read/*prin*
	done
}

function make_SAF {
	local saf_name=$1.saf
	local fasta=$2
	local source=$3
	local feature=$4		

	> tmp
	paste <(awk -F'>' '{print $2}' $fasta) tmp | column -s $'\t' -t > tmp
	sed -i "s/$/\t$(echo $source)/" tmp
	sed -i "s/$/\t$(echo $feature)/" tmp
	sed -i "s/$/\t$(echo 1)/" tmp
	paste tmp <(egrep "^[^>]" $fastalocal path_data_fastq=/mnt/data/troyalty_data/sequencing/fastq_files/ | awk '{print length}') | column -s $'\t' -t > $saf_name
	rm tmp
	sed -i "s/$/\t$(echo .)/" $saf_name
 
}

function single_line {
        local file_path=$1
        local new_file=$2
        local line_flag=
        local merge_flag=0
#	local count=/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $path_hold$read_tailer_single > $path_hold$read_tailer_single$read_tailer_short
                cat $path_hold$read_tailer_short0

        > $new_file

        while read line; do
                echo $line | grep -q ">"
                if [[ $? -eq 0 ]]; then
                        echo $line >> $new_file
                        merge_flag=0
#			let count=$count+1
                elif [[ $merge_flag -eq 0 ]]; then
                        echo $line >> $new_file
                        merge_flag=1
#			let count=$count+1
                else
#			sed -i "${count}N;s/\n//" $file_path
                   	echo "$(cat $new_file)$line" > $new_file
                fi
        done < $file_path
}


function rename_contig {
	local read_array=($@)
	local read_header=SRR
	local read_tailer=_short.fasta
	local read_rename=_rename.fasta
	local path_data_megahit=/mnt/data/troyalty_data/sequencing/megahit_contig/
	local path_hold=
	local file_name=
	local count1=1
	local count2=1
	local check=

	for read_num in $read_array; do
		file_name=$path_data_megahit$read_header$read_num/$read_header$read_num$read_tailer
		path_hold=$path_data_megahit$read_header$read_num/$read_header$read_num$read_rename 
		> $path_hold
			while read line; do
				check=$(bc <<< "$count1%2")
				if [[ $check -eq 1 ]]; then
					contig_tail=
					echo ">$read_header$read_num/$count2" >> $path_hold
					let count2="count2 + 1"
				else
					echo $line >> $path_hold
				fi
				let count1="count1 + 1"
		
			done < $file_name
		count=1	
	done
}

#Downloads specified reads from the configuration file
function download_reads {
	local read_array=$@
	local read_header=SRR
	local path_data=/mnt/data/troyalty_data/sequencing/fastq_files/

	for read in $read_array; do
		rm -fr $path_data$read_header$read
		mkdir $path_data$read_header$read
		echo Downloading $read_header$read
		parallel-fastq-dump -s $read_header$read -t 20 -O $path_data$read_header$read --split-files
#		fastq-dump --split-files $read_header$read --outdir $path_data$read_header$read
#		rm -fr $path_data/$read_header$read
#		mv $read_header$read* $path_data/$read_header$read
	done 
}



#Remove adapters from fastq files using Trimmomatic
function apply_trimm {
	local array=($@)
	local num_thread=$1
	local read_header=SRR
	local read_tailer_forward_frac=_fraction_1.fastq
	local read_tailer_reverse_frac=_fraction_2.fastq
	local read_tailer_forward=_1.fastq
        local read_tailer_reverse=_2.fastq
	local read_paired=_paired_output
	local read_unpaired=_unpaired_output
	local path_data_fastq=/mnt/data/troyalty_data/sequencing/fastq_files/
	local path_data_fastq_sub=
	local path_data_fastq_trim=/mnt/data/troyalty_data/sequencing/trimmed_fastq_files/
	local path_trimm=~/Software/Trimmomatic-0.36/trimmomatic-0.36.jar
	local path_data_trimm_fastq=/mnt/data/troyalty_data/sequencing/trimmed_fastq_files/
	local path_adap=~/Software/Trimmomatic-0.36/adapters/TruSeq2-PE.fa
	local current_file_forward=
	local current_file_reverse=
	local output_forward_paired=
	local output_reverse_paired=
	local output_forward_unpaired=
	local output_reverse_unpaired=
	local i=1
#	local max_it=$(echo "$# - 1" | bc -l)

	while [[ $i -le $#-1 ]]; do
		read=${array[$i]}
		path_data_fastq_sub=$path_data_fastq$read_header$read
		current_file_forward=$path_data_fastq_sub/$read_header$read$read_tailer_forward_frac
		current_file_reverse=$path_data_fastq_sub/$read_header$read$read_tailer_reverse_frac
		output_forward_paired=$path_data_fastq_trim/$read_header$read/$read_header$read$read_paired$read_tailer_forward
		output_reverse_paired=$path_data_fastq_trim/$read_header$read/$read_header$read$read_paired$read_tailer_reverse
		output_forward_unpaired=$path_data_fastq_trim/$read_header$read/$read_header$read$read_unpaired$read_tailer_forward
		output_reverse_unpaired=$path_data_fastq_trim/$read_header$read/$read_header$read$read_unpaired$read_tailer_reverse

		rm -fr $path_data_fastq_trim/$read_header$read
		mkdir $path_data_fastq_trim/$read_header$read

		java -jar $path_trimm PE -threads $num_thread $current_file_forward $current_file_reverse $output_forward_paired $output_forward_unpaired $output_reverse_paired $output_reverse_unpaired ILLUMINACLIP:$path_adap:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36	
		let i=$i+1
	done
}

function apply_mega {
	local read_array=($@)
	local read_header=SRR
	local read_tailer=.fastq
	local read_trim_tailer_1=_paired_output_1
	local read_trim_tailer_2=_paired_output_2
	local path_megahit_out=/mnt/data/troyalty_data/sequencing/megahit_contig/
	local path_data_trimm=/mnt/data/troyalty_data/sequencing/trimmed_fastq_files/
	

	for read in $read_array; do
		rm -fr $path_megahit_out$read_header$read
		/home/troyalty/Software/megahit-master/megahit -1 $path_data_trimm$read_header$read/$read_header$read$read_trim_tailer_1$read_tailer -2 $path_data_trimm$read_header$read/$read_header$read$read_trim_tailer_2$read_tailer -o $path_megahit_out$read_header$read --presets meta-large -t 32
#--k-min 45 --k-max 85 --k-step 10
#		megahit -1 $path_data_trimm$read_header$read/$read_header$read$read_trim_tailer_1$read_tailer -2 $path_data_trimm$read_header$read/$read_header$read$read_trim_tailer_2$read_tailer -o $path_megahit_out$read_header$read --preset meta-sensitive
	done
}


function remove_overlap {
	local read_array=($@)
	local read_header=SRR
	local read_tailer=_short.fasta
	local read_tailer_seq=_seq.fasta
	local read_tailer_hold=_hold
	local read_tailer_hold2=_hold2
	local path_data_megahit=/mnt/data/troyalty_data/sequencing/megahit_contig/
	local file_contig=final.contigs.fa	

	for read in $read_array; do
		seqmagick convert --min-length $contig_min $path_data_megahit$read_header$read/$file_contig $path_data_megahit$read_header$read/$read_header$read$read_tailer_seq
		echo $path_data_megahit$read_header$read/$file_contig
		cd-hit-est -i $path_data_megahit$read_header$read/$read_header$read$read_tailer_seq -o $path_data_megahit$read_header$read/$read_header$read$read_tailer_hold -T 90 -M 500000 -c 0.99 -n 10
		awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $path_data_megahit$read_header$read/$read_header$read$read_tailer_hold > $path_data_megahit$read_header$read/$read_header$read$read_tailer_hold2
		awk -v var=">$read_header$read" 'BEGIN { cntr = 0} /^>/ { cntr++  ; print var , cntr } !/^>/ { print $0 }' $path_data_megahit$read_header$read/$read_header$read$read_tailer_hold2 > $path_data_megahit$read_header$read/$read_header$read$read_tailer
		sed -i "s/\s/_/" $path_data_megahit$read_header$read/$read_header$read$read_tailer
		rm $path_data_megahit$read_header$read/*hold*
#		rename_contig ${read_array[@]}
	done
}


function amos_assemble {
	local read_array=($@)
	local read_header=SRR
	local read_tailer_megahit=_short.fasta
	local read_tailer_single=.singletons.seq
	local read_tailer_megahit_afg=_short.fasta.afg
	local read_tailer_hold=_hold
	local read_tailer_second=_SECONDARY_contigs.fasta
	local read_tailer_fasta=.fasta
	local read_tailer_combined=.comb
	local path_amos=~/Software/amos-3.1.0/bin/
	local path_data_megahit=/mnt/data/troyalty_data/sequencing/megahit_contig/	
	local path_data_amos=/mnt/data/troyalty_data/sequencing/amos/
	local path_hold=
	local path_afg=

	for read in $read_array; do
		path_hold=$path_data_megahit$read_header$read/$read_header$read$read_tailer_megahit
		if [[ $multiple_sites -eq 1 ]]; then
			path_afg=$path_data_megahit$read_header$read/$read_header$read$read_tailer_megahit_afg
			rm -fr $path_data_amos$read_header$read
			mkdir $path_data_amos$read_header$read
			toAmos -s $path_hold -o $path_afg
			minimus2 $path_hold -D OVERLAP=100 MINID=95
			awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $path_hold$read_tailer_fasta  > $path_hold$read_tailer_combined
			awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $path_hold$read_tailer_single > $path_hold$read_tailer_single$read_tailer_hold
			sed -i "s/>.*/&_PRIMARY/" $path_hold$read_tailer_single$read_tailer_hold 
			cat $path_hold$read_tailer_combined $path_hold$read_tailer_single$read_tailer_hold > $path_data_amos$read_header$read/$read_header$read$read_tailer_second
			rm $path_data_megahit$read_header$read/*hold*
##			rm $path_hold$read_tailer_single$read_tailer_short
			mv $path_hold.* $path_data_amos$read_header$read/		
		else
			rm -fr $path_data_amos$read_header$read
			mkdir $path_data_amos$read_header$read
			cp $path_hold $path_data_amos$read_header$read/$read_header$read$read_tailer_second
		fi
	done
}


function apply_bowtie2 {
	local read_array=($@)
	local read_header=SRR
	local read_tailer_amos=_SECONDARY_contigs.fasta
	local read_tailer_bowtie=_SECONDARY_contigs.bt_index
	local read_tailer_id=.id
	local read_tailer_sam=.sam
#	local reader_tailer_saf=.saf
	local read_tailer_trim_forward=_paired_output_1.fastq
	local read_tailer_trim_backward=_paired_output_2.fastq
	local path_data_amos=/mnt/data/troyalty_data/sequencing/amos/
	local path_data_bowtie=/mnt/data/troyalty_data/sequencing/bowtie_assemblies/
	local path_data_trimm=/mnt/data/troyalty_data/sequencing/trimmed_fastq_files/
	local threads=$1
	local i=1
	
	while [[ $i -le $#-1 ]]; do
		read=${read_array[$i]}
		rm -fr $path_data_bowtie$read_header$read 
		mkdir $path_data_bowtie$read_header$read
		bowtie2-build $path_data_amos$read_header$read/$read_header$read$read_tailer_amos $path_data_bowtie$read_header$read/$read_header$read$read_tailer_bowtie
		bowtie2 -q -1 $path_data_trimm$read_header$read/$read_header$read$read_tailer_trim_forward -2 $path_data_trimm$read_header$read/$read_header$read$read_tailer_trim_backward -S $path_data_bowtie$read_header$read/$read_header$read$read_tailer_sam -x $path_data_bowtie$read_header$read/$read_header$read$read_tailer_bowtie --no-unal -p $threads
		let i=$i+1
	done
}


function apply_binsanity {
	local read_array=($@)
	local read_header=SRR
	local read_tailer_amos=_SECONDARY_contigs.fasta	
	local read_tailer_bam=.bam
	local read_tailer_sam=.sam
	local read_tailer_id=.ids
	local read_tailer_fna=*.fna
	local read_tailer_cov=.cov.x100.lognorm
	local path_data_amos=/mnt/data/troyalty_data/sequencing/amos/
	local path_data_bowtie=/mnt/data/troyalty_data/sequencing/bowtie_assemblies/
	local path_data_binsanity=/mnt/data/troyalty_data/sequencing/binsanity/
	local path_data_results=BINSANITY-RESULTS/
	local read_cutoff=$1
	local i=1
	local count2=
	local trip_header=trip_
	local tail=
	local outfile=taxonomic_summary	
	local tmp_dir=

	while [[ $i -le $#-1 ]]; do	
		read=${read_array[$i]}
		rm -fr $path_data_binsanity$read_header$read
		mkdir $path_data_binsanity$read_header$read
		cd $path_data_binsanity$read_header$read
		samtools view -bS $path_data_bowtie$read_header$read/$read_header$read$read_tailer_sam > $path_data_binsanity$read_header$read/$read_header$read$read_tailer_bam
		get-ids -f $path_data_amos$read_header$read -l $path_data_amos$read_header$read/$read_header$read$read_tailer_amos -o $path_data_binsanity$read_header$read/$read_header$read$read_tailer_id  -x $read_cutoff
		Binsanity-profile -i $path_data_amos$read_header$read/$read_header$read$read_tailer_amos -s $path_data_binsanity$read_header$read/ --ids $path_data_binsanity$read_header$read/$read_header$read$read_tailer_id -c $path_data_binsanity$read_header$read/$read_header$read
		cp $path_data_amos$read_header$read/$read_header$read$read_tailer_amos $path_data_binsanity$read_header$read

		cd $path_data_binsanity$read_header$read/
		custom_binsanity_wf $outfile 32 $path_data_binsanity$read_header$read/$read_header$read$read_tailer_amos
		
		nucleotide_to_protein bins
		parallel_rps_blast bins
		
		tmp_dir=$site_dir/$fraction
		if [[ -d "$site_dir" ]]; then
			if [[ -d "$tmp_dir" ]]; then
				count2=$(ls $site_dir/$fraction | wc -l | cut -d ' ' -f 1)
				count2=$(echo "$count2+1" | bc -l)
				mkdir $site_dir/$fraction/$trip_header$count2
			else
				count2=1	
				mkdir $site_dir/$fraction
				mkdir $site_dir/$fraction/$trip_header$count2
			fi
		else 
			count2=1
			mkdir $site_dir
			mkdir $site_dir/$fraction
                        mkdir $site_dir/$fraction/$trip_header$count2
		fi
		mv bins $outfile *cov.x100.lognorm *SECONDARY* $tmp_dir/$trip_header$count2

		let i=$i+1
	done

}

##########
#Script
##########

org_dir=$(pwd)
cd ~/Software/assembly_pipeline

source configuration_files/$config_file
source ~/Software/assembly_pipeline/checkm_wf_taxa.sh
source ~/Software/assembly_pipeline/genomic_functions.sh
source ~/Software/assembly_pipeline/compile_csv.sh

if [[ $download -eq 1 ]]; then
	download_reads ${reads[@]}
fi

if [[ $dedupe -eq 1 ]]; then
	remove_fastq_duplicates ${reads[@]}
fi

for iteration in $(seq 1 $num_it); do

        if [[ $sub -eq 1 ]]; then
                randomize_contigs ${reads[@]}
        fi

	if [[ $trimm -eq 1 ]]; then
		apply_trimm $threads ${reads[@]}
	fi

	if [[ $mega -eq 1 ]]; then
		apply_mega ${reads[@]}
	fi

	if [[ $overlap -eq 1 ]]; then
		remove_overlap ${reads[@]}
	fi

	if [[ $amos -eq 1 ]]; then
		amos_assemble ${reads[@]}	
	fi

	if [[ $bowtie -eq 1 ]]; then
		apply_bowtie2 $threads ${reads[@]}	
	fi

	if [[ $binsanity -eq 1 ]]; then
		apply_binsanity $contig_min ${reads[@]}
	fi

	cd $org_dir
done
