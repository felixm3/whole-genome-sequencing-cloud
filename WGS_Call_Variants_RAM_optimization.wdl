import "https://api.firecloud.org/ga4gh/v1/tools/gatk:alignment/versions/4/plain-WDL/descriptor" as Alignment
import "https://api.firecloud.org/ga4gh/v1/tools/gatk:split-large-readgroup/versions/4/plain-WDL/descriptor" as SplitRG
import "https://api.firecloud.org/ga4gh/v1/tools/gatk:quality-control/versions/2/plain-WDL/descriptor" as QC
import "https://api.firecloud.org/ga4gh/v1/tools/gatk:bam-processing/versions/4/plain-WDL/descriptor" as Processing
import "https://api.firecloud.org/ga4gh/v1/tools/gatk:germline-variant-discovery/versions/3/plain-WDL/descriptor" as Calling
import "https://api.firecloud.org/ga4gh/v1/tools/gatk:utilities/versions/2/plain-WDL/descriptor" as Utils
import "https://raw.githubusercontent.com/gatk-workflows/gatk4-cnn-variant-filter/1.2.0/tasks/cnn_variant_common_tasks.wdl" as CNNTasks

# WORKFLOW DEFINITION

workflow WGS_Call_Variants {

		String SRA_ID

		String bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta"

		File ref_fasta
		File ref_fasta_index
		File ref_dict
		File ref_alt
		File ref_bwt
		File ref_sa
		File ref_amb
		File ref_ann
		File ref_pac

        File contamination_sites_ud
        File contamination_sites_bed
        File contamination_sites_mu

        File dbSNP_vcf
        File dbSNP_vcf_index
        Array[File] known_indels_sites_VCFs
        Array[File] known_indels_sites_indices

        File wgs_calling_interval_list

        String gatk_docker

        Array[File] resources            # List of VCF file names of resources of known SNPs and INDELs, (e.g. mills, gnomAD)
        Array[File] resources_index      # List of VCF file indices of resources




	# download data from SRA and convert to uBAM
	call Fetch_SRA_to_uBAM {
        input:
            SRA_ID = SRA_ID
    }

	# Get the version of BWA to include in the PG record in the header of the BAM produced
	# by MergeBamAlignment.
	call Alignment.GetBwaVersion

	Int compression_level = 2
	Int preemptible_tries = 3


    # Split bam into multiple smaller bams,
    # map reads to reference and recombine into one bam
    call SplitRG.split_large_readgroup as SplitRG {
        input:
          input_bam = Fetch_SRA_to_uBAM.unmapped_bam,
          bwa_commandline = bwa_commandline,
          bwa_version = GetBwaVersion.version,
          output_bam_basename = basename(Fetch_SRA_to_uBAM.unmapped_bam, ".bam") + ".aligned.unsorted",
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          ref_alt = ref_alt,
          ref_amb = ref_amb,
          ref_ann = ref_ann,
          ref_bwt = ref_bwt,
          ref_pac = ref_pac,
          ref_sa = ref_sa,
          compression_level = compression_level,
          preemptible_tries = preemptible_tries
    }

    File output_aligned_bam = SplitRG.aligned_bam

    Float mapped_bam_size = size(output_aligned_bam, "GB")

    # Aggregate aligned+merged flowcell BAM files and mark duplicates
    # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
    # to avoid having to spend time just merging BAM files.
    call Processing.MarkDuplicates as MarkDuplicates {
        input:
          input_bams = [output_aligned_bam],
          output_bam_basename = SRA_ID + ".aligned.unsorted.duplicates_marked",
          metrics_filename = SRA_ID + ".duplicate_metrics",
          total_input_size = mapped_bam_size,
          compression_level = compression_level,
          preemptible_tries = preemptible_tries
    }

    # Sort aggregated+deduped BAM file
    call Processing.SortSam as SortSampleBam {
        input:
          input_bam = MarkDuplicates.output_bam,
          output_bam_basename = SRA_ID + ".aligned.duplicate_marked.sorted",
          compression_level = compression_level,
          preemptible_tries = preemptible_tries
    }

    # Create list of sequences for scatter-gather parallelization
    call Utils.CreateSequenceGroupingTSV as CreateSequenceGroupingTSV {
        input:
          ref_dict = ref_dict,
          preemptible_tries = preemptible_tries
    }

    # Estimate level of cross-sample contamination
    call Processing.CheckContamination as CheckContamination {
        input:
          input_bam = SortSampleBam.output_bam,
          input_bam_index = SortSampleBam.output_bam_index,
          contamination_sites_ud = contamination_sites_ud,
          contamination_sites_bed = contamination_sites_bed,
          contamination_sites_mu = contamination_sites_mu,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          output_prefix = SRA_ID + ".preBqsr",
          preemptible_tries = preemptible_tries,
          contamination_underestimation_factor = 0.75
    }

    # BQSR
    # We need disk to localize the sharded input and output due to the scatter for BQSR.
    # If we take the number we are scattering by and reduce by 3 we will have enough disk space
    # to account for the fact that the data is not split evenly.
    Int num_of_bqsr_scatters = length(CreateSequenceGroupingTSV.sequence_grouping)
    Int potential_bqsr_divisor = num_of_bqsr_scatters - 10
    Int bqsr_divisor = if potential_bqsr_divisor > 1 then potential_bqsr_divisor else 1

    # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
    scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
        # Generate the recalibration model by interval
        call Processing.BaseRecalibrator as BaseRecalibrator {
          input:
            input_bam = SortSampleBam.output_bam,
            recalibration_report_filename = SRA_ID + ".recal_data.csv",
            sequence_group_interval = subgroup,
            dbSNP_vcf = dbSNP_vcf,
            dbSNP_vcf_index = dbSNP_vcf_index,
            known_indels_sites_VCFs = known_indels_sites_VCFs,
            known_indels_sites_indices = known_indels_sites_indices,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            bqsr_scatter = bqsr_divisor,
            preemptible_tries = preemptible_tries
        }
    }

    # Merge the recalibration reports resulting from by-interval recalibration
    # The reports are always the same size
    call Processing.GatherBqsrReports as GatherBqsrReports {
    input:
      input_bqsr_reports = BaseRecalibrator.recalibration_report,
      output_report_filename = SRA_ID + ".recal_data.csv",
      preemptible_tries = preemptible_tries
    }

    scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {
        # Apply the recalibration model by interval
        call Processing.ApplyBQSR as ApplyBQSR {
          input:
            input_bam = SortSampleBam.output_bam,
            output_bam_basename = SRA_ID,
            recalibration_report = GatherBqsrReports.output_bqsr_report,
            sequence_group_interval = subgroup,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            bqsr_scatter = bqsr_divisor,
            compression_level = compression_level,
            preemptible_tries = preemptible_tries
        }
    }

    Float agg_bam_size = size(SortSampleBam.output_bam, "GB")

    # Merge the recalibrated BAM files resulting from by-interval recalibration
    call Processing.GatherSortedBamFiles as GatherBamFiles {
        input:
          input_bams = ApplyBQSR.recalibrated_bam,
          output_bam_basename = SRA_ID,
          total_input_size = agg_bam_size,
          compression_level = compression_level,
          preemptible_tries = preemptible_tries
    }

    #BQSR bins the qualities which makes a significantly smaller bam
    Float binned_qual_bam_size = size(GatherBamFiles.output_bam, "GB") 

    # ValidateSamFile runs out of memory in mate validation on crazy edge case data, so we want to skip the mate validation
    # in those cases.  These values set the thresholds for what is considered outside the normal realm of "reasonable" data.
    Float max_duplication_in_reasonable_sample = 0.30
    Float max_chimerism_in_reasonable_sample = 0.15

    # Convert the final merged recalibrated BAM file to CRAM format
    call Utils.ConvertToCram as ConvertToCram {
        input:
          input_bam = GatherBamFiles.output_bam,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          output_basename = SRA_ID,
          preemptible_tries = preemptible_tries
    }

    # Validate the CRAM file
    call QC.ValidateSamFile as ValidateCram {
        input:
          input_bam = ConvertToCram.output_cram,
          input_bam_index = ConvertToCram.output_cram_index,
          report_filename = SRA_ID + ".cram.validation_report",
          ref_dict = ref_dict,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ignore = ["MISSING_TAG_NM"],
          max_output = 1000000000,
          is_outlier_data = false,
          preemptible_tries = preemptible_tries
    }

    Int scatter_count = 50
    Int additional_disk = 20
    Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")

    # Break the calling interval_list into sub-intervals
    # Perform variant calling on the sub-intervals, and then gather the results
    call Utils.ScatterIntervalList as ScatterIntervalList {
        input:
          interval_list = wgs_calling_interval_list,
          scatter_count = scatter_count,
          break_bands_at_multiples_of = 1000000
    }

    # We need disk to localize the sharded input and output due to the scatter for HaplotypeCaller.
    # If we take the number we are scattering by and reduce by 20 we will have enough disk space
    # to account for the fact that the data is quite uneven across the shards.
    Int potential_hc_divisor = ScatterIntervalList.interval_count - 20
    Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1

    # Call variants in parallel over WGS calling intervals
    scatter (index in range(ScatterIntervalList.interval_count)) {
        call Calling.HaplotypeCaller_GATK4_VCF as HaplotypeCaller4 {
            input:
                contamination = CheckContamination.contamination,
                input_bam = GatherBamFiles.output_bam,
                interval_list = ScatterIntervalList.out[index],
                make_gvcf = false,
                vcf_basename = SRA_ID,
                ref_dict = ref_dict,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                hc_scatter = hc_divisor,
                preemptible_tries = preemptible_tries
        }

        call CNNTasks.CNNScoreVariants {
            input:
                input_vcf = HaplotypeCaller4.output_vcf,
                input_vcf_index = HaplotypeCaller4.output_vcf_index,
                reference_fasta = ref_fasta,
                reference_dict = ref_dict,
                reference_fasta_index = ref_fasta_index,               
                output_prefix = SRA_ID,
                interval_list = ScatterIntervalList.out[index],
                gatk_docker = gatk_docker,
                preemptible_attempts = preemptible_tries,
                mem_gb = 7,
                disk_space_gb = round((binned_qual_bam_size/scatter_count) + ref_size + additional_disk)
        }

     File merge_input = HaplotypeCaller4.output_vcf
     File merge_input_index = HaplotypeCaller4.output_vcf_index
    }


    # Combine by-interval raw VCFs into a single sample VCF file
    call Calling.MergeVCFs as MergeVCFs {
        input:
              input_vcfs = merge_input,
              input_vcfs_indexes = merge_input_index,
              output_vcf_name = SRA_ID + ".raw.vcf.gz",
              preemptible_tries = preemptible_tries
    }

    # Combine by-interval CNN VCFs into a single sample VCF file
    call CNNTasks.MergeVCFs as MergeVCFsCNN {
        input: 
            input_vcfs = CNNScoreVariants.cnn_annotated_vcf,
            output_prefix = SRA_ID,
            preemptible_attempts = preemptible_tries,
            gatk_docker = gatk_docker,
            disk_space_gb = additional_disk
    }

    call CNNTasks.FilterVariantTranches {
        input:
            input_vcf = MergeVCFsCNN.merged_vcf,
            input_vcf_index = MergeVCFsCNN.merged_vcf_index,
            resources = resources,
            resources_index = resources_index,
            output_prefix = SRA_ID,
            snp_tranches = " --snp-tranche 99.9 ",
            indel_tranches = " --indel-tranche 99.5 ",
            info_key = "CNN_1D",
            preemptible_attempts = preemptible_tries,
            gatk_docker = gatk_docker,
            disk_space_gb = additional_disk
    }

    # Outputs that will be retained when execution is complete
    output {
        File selfSM = CheckContamination.selfSM
        Float contamination = CheckContamination.contamination

        File duplicate_metrics = MarkDuplicates.duplicate_metrics
        File output_bqsr_reports = GatherBqsrReports.output_bqsr_report

        File output_cram = ConvertToCram.output_cram
        File output_cram_index = ConvertToCram.output_cram_index
        File output_cram_md5 = ConvertToCram.output_cram_md5

        File validate_cram_file_report = ValidateCram.report

        File output_vcf_raw = MergeVCFs.output_vcf
        File output_vcf_raw_index = MergeVCFs.output_vcf_index

        File output_vcf_CNN_filtered = FilterVariantTranches.cnn_filtered_vcf
        File output_vcf_CNN_filtered_index = FilterVariantTranches.cnn_filtered_vcf_index 
    }



}

task Fetch_SRA_to_uBAM {

    String  SRA_ID
    Int?    machine_mem_gb
    String  docker = "quay.io/broadinstitute/ncbi-tools:2.10.7.10"

    Int disk_size = 750
    meta {
        description: "This searches NCBI SRA for accessions using the Entrez interface, collects associated metadata, and returns read sets as unaligned BAM files with metadata loaded in. directly."
    }
    command {
        set -e
        # fetch SRA metadata on this record
        esearch -db sra -q ${SRA_ID} | efetch -mode json -json > SRA.json
        cp SRA.json ${SRA_ID}.json

        # pull reads from SRA and make a fully annotated BAM -- must succeed
        PLATFORM=$(jq -r '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.EXPERIMENT.PLATFORM | keys[] as $k | "\($k)"' SRA.json)
        SAMPLE=$(jq -r '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.SAMPLE.IDENTIFIERS.EXTERNAL_ID|select(.namespace == "BioSample")|.content' SRA.json)
        LIBRARY=$(jq -r .EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.EXPERIMENT.alias SRA.json)

        sam-dump --unaligned --header ${SRA_ID} \
            | samtools view -bhS - \
            > temp.bam
        
        picard AddOrReplaceReadGroups \
            I=temp.bam \
            O=${SRA_ID}.bam \
            RGLB="$LIBRARY" \
            RGSM="$SAMPLE" \
            RGPL="$PLATFORM" \
            RGPU="$LIBRARY" \
            SORT_ORDER=queryname \
            MAX_RECORDS_IN_RAM=2000000 \
            -Xmx64G \
            VALIDATION_STRINGENCY=SILENT
        rm temp.bam

        samtools view -H ${SRA_ID}.bam

        # emit numeric WDL outputs
        echo $PLATFORM > OUT_PLATFORM
        echo $SAMPLE > OUT_BIOSAMPLE
        echo $LIBRARY > OUT_LIBRARY
        samtools view -c ${SRA_ID}.bam | tee OUT_NUM_READS
    }

    output {
        File    unmapped_bam              = "${SRA_ID}.bam"
        Int     num_reads                 = read_int("OUT_NUM_READS")
        String  sequencing_platform       = read_string("OUT_PLATFORM")
        String  biosample_accession       = read_string("OUT_BIOSAMPLE")
        String  library_id                = read_string("OUT_LIBRARY")
        File    sra_metadata              = "${SRA_ID}.json"

    }

    runtime {
        cpu:     16
        memory:  64 + " GB"
        disks:   "local-disk " + disk_size + " LOCAL"
        disk:    disk_size + " GB" # TES
        dx_instance_type: "mem2_ssd1_v2_x2"
        docker:  docker
        #maxRetries: 2
    }
}
