'''

master script for running pipeline to call variants on WES tumor data
mpnst samples
written by ben demaree 10.21.2018

file requirements:
-paired fastq files for each tumor sample (multiple lanes)
-human reference genome (bt2 index, fasta, dict, fai)

software requirements:
-bowtie2
-samtools
-picard
-gatk

'''

import os
import subprocess
from slackclient import SlackClient

class TumorSample(object):
    # class for storing metadata for each tumor sample

    def __init__(self, fastq_r1, fastq_r2, sample_name):
        # initialize object with list of fastq files and sample name

        self.fastq_r1 = fastq_r1                # list of fastq files for this sample (R1)
        self.fastq_r2 = fastq_r2                # list of fastq files for this sample (R2)
        self.sample_name = sample_name          # base name of sample

    def generate_filenames(self):
        # generate all filenames for pipeline

        self.sample_folder = out_dir + self.sample_name + '/'                           # folder for this sample
        self.concat_R1 = self.sample_folder + self.sample_name + '_R1.fastq.gz'         # R1 concatenated file
        self.concat_R2 = self.sample_folder + self.sample_name + '_R2.fastq.gz'         # R2 concatenated file
        self.bowtie_stats = self.sample_folder + self.sample_name + '_bowtiestats.txt'  # bowtiestats file
        self.bamfile = self.sample_folder + self.sample_name + '.bam'                   # bam file
        self.bamfile_dedup = add_ext(self.bamfile, 'dedup')                             # deduplicated bam file
        self.dedup_metrics = self.sample_folder + self.sample_name + '.dedup.metrics'   # deduplication metrics file
        self.bamfile_index = self.bamfile_dedup[:-4] + '.bai'                           # bam file index
        self.vcffile = self.sample_folder + self.sample_name + '.vcf.gz'                # vcf file
        self.realigned_bam = self.sample_folder + self.sample_name + '.realigned.bam'   # gatk-realigned bam file

    def create_sample_folder(self):
        # create data folder for a sample

        os.mkdir(self.sample_folder)

    def concatenate_fastq(self):
        # concatenate fastqs from multiple lanes for a given sample

        concat_cmd = 'cat %s > %s; cat %s > %s' \
                     % (' '.join(self.fastq_r1),
                        self.concat_R1,
                        ' '.join(self.fastq_r2),
                        self.concat_R2)

        process = subprocess.Popen(concat_cmd, shell=True)

        return process

    def align_sample(self):
        # align the panel to the bowtie2 human index

        align_cmd = 'cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 40 --quiet --interleaved %s %s' \
                    ' | (/usr/local/bin/bowtie2-2.3.4.1-linux-x86_64/bowtie2 -x %s --interleaved -' \
                    ' --rg-id %s --rg SM:%s --rg PL:ILLUMINA --rg CN:UCSF) 2>%s' \
                    ' | samtools view -b' \
                    ' | samtools sort -o %s -m 100M' \
                     % (self.concat_R1,
                        self.concat_R2,
                        bt2_ref,
                        self.sample_name,
                        self.sample_name,
                        self.bowtie_stats,
                        self.bamfile)

        process = subprocess.Popen(align_cmd, shell=True)

        return process

    def mark_duplicates(self):
        # mark duplicates in all bam files using picard

        md_cmd = 'picard MarkDuplicates I=%s O=%s M=%s' \
                % (self.bamfile,
                   self.bamfile_dedup,
                   self.dedup_metrics)

        process = subprocess.Popen(md_cmd, shell=True)

        return process

    def index_bam(self):
        # index all bam files using samtools

        index_cmd = 'samtools index %s %s' \
                    % (self.bamfile_dedup,
                       self.bamfile_index)

        process = subprocess.Popen(index_cmd, shell=True)

        return process

    def call_variants(self):
        # call variants using gatk

        variants_cmd = 'gatk HaplotypeCaller -R %s -I %s -O %s -ERC GVCF -L %s -D %s -bamout %s' \
                       ' --interval-padding 100 --native-pair-hmm-threads 10' \
                    % (human_fasta,
                       self.bamfile_dedup,
                       self.vcffile,
                       exome_bed,
                       dbsnp_file,
                       self.realigned_bam)

        process = subprocess.Popen(variants_cmd, shell=True)

        return process

def file_summary(samples, filename='file_summary.txt', to_file=True):
    # displays sample files and class variables

    if to_file:
        f = open(filename, 'w')

    for sample in samples:

        s = vars(sample)

        for item in s:
            print item,': ', s[item]
        print '\n'

        if to_file:
            for item in s:
                if type(s[item]) is list:
                    f.write(item + ': ' + ', '.join(s[item]) + '\n')

                else:
                    f.write(item + ': ' + s[item] + '\n')
            f.write('\n')

    if to_file:
        f.close()

    print '%d samples identified.\n' % len(samples)

def get_fastq_filenames(fastq_dir, paired=True):
    # gets fastq filenames in a given directory and runs some simple checks
    # assumes files are compressed with .gz extensions

    R1_files = []
    R2_files = []

    for dir in fastq_dir:

        for file in os.listdir(dir):

            if file.endswith('.fastq.gz'):
                R = file.split('_')[-2]  # R1 or R2

                if R == 'R1':
                    R1_files += [dir + file]

                elif R == 'R2' and paired:
                    R2_files += [dir + file]

                else:
                    print 'Unexpected filename structure. Exiting...'
                    raise SystemExit

    if len(R1_files) != len(R2_files) and paired:
        print 'Unpaired FASTQ files exist! Check input files.'
        raise SystemExit

    R1_files.sort()
    R2_files.sort()

    return R1_files, R2_files

def partition_fastq(R1_files, R2_files):
    # assigns fastq files to a TumorSample object

    sample_dict = {}

    for i in range(len(R1_files)):

        sample_name = R1_files[i].split('/')[-1].split('_')[0]

        try:
            sample_dict[sample_name]['R1'].append(R1_files[i])
            sample_dict[sample_name]['R2'].append(R2_files[i])

        except KeyError:
            sample_dict[sample_name] = {}
            sample_dict[sample_name]['R1'] = []
            sample_dict[sample_name]['R2'] = []

            sample_dict[sample_name]['R1'].append(R1_files[i])
            sample_dict[sample_name]['R2'].append(R2_files[i])

    # store in TumorSample objects
    samples = []

    for s in sample_dict:
        samples.append(TumorSample(sample_dict[s]['R1'], sample_dict[s]['R2'], s))

    return samples

def add_ext(filename, string):
    # adds a string before the extension in a filename (do not add period to input string)

    return '.'.join(['.'.join(filename.split('.')[:-1]), str(string), filename.split('.')[-1]])

def wait(processes):
    # waits for processes to finish

    return [process.communicate() for process in processes]

def slack_message(message):
    # for posting a notification to the server-alerts slack channel

    channel = 'server-alerts'
    token = 'xoxp-7171342752-7171794564-486340412737-91fd92781cde6307b077f30f9ea1b700'
    sc = SlackClient(token)

    sc.api_call('chat.postMessage', channel=channel,
                text=message, username='pipelines',
                icon_emoji=':adam:')

if __name__ == "__main__":

    # global experiment variables - modify for each batch of samples

    # output directory base
    out_dir = '/drive2/mpnst_out/'

    # summary file path
    summary_file = out_dir + 'mpnst_sample_summary.txt'

    # input fastq file directories
    fastq_dir = ['/drive2/hvasu/DR11_Fastq/',
                 '/drive2/hvasu/DR12_Fastq/',
                 '/drive2/hvasu/DR13_Fastq/',
                 '/drive2/hvasu/DR14_Fastq/']

    # bowtie2 index location
    bt2_ref = '/drive2/igenomes/hg19/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome'

    # human reference genome fasta file path
    human_fasta = '/drive2/igenomes/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa'

    # exome coordinates file path
    exome_bed = '/drive2/hvasu/SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_hg19_capture_targets.bed'

    # dbsnp file path for variant calling
    dbsnp_file = '/drive2/hvasu/common_all_20180423.vcf.gz'

    print '''

Beginning tumor WES variant calling pipeline...

####################################################################################
# Step 0: identify samples
####################################################################################
    '''

    # get all fastq filenames
    R1_files, R2_files = get_fastq_filenames(fastq_dir)

    # partition fastq files and initialize TumorSample objects
    samples = partition_fastq(R1_files, R2_files)

    # generate filenames for each sample
    for sample in samples:
        sample.generate_filenames()

    # make folder for each sample
    for sample in samples:
        try:
            sample.create_sample_folder()
        except OSError:
            print 'Folder %s already exists. Please move or rename the existing folder!' % \
                  (sample.sample_folder)
            raise SystemExit

    # display and write sample summary file
    file_summary(samples, summary_file, to_file=True)

    print '''
####################################################################################
# Step 1: concatenate fastq files
####################################################################################
'''

    # concatenate fastq files for each sample
    concatenate = [sample.concatenate_fastq() for sample in samples]
    # wait for all processes to finish before continuing
    wait(concatenate)

    print '''
####################################################################################
# Step 2: pre-process files (cut adapters, align, convert to bam, sort)
####################################################################################
'''

    # perform pre-processing
    pre_process = [sample.align_sample() for sample in samples]
    # wait for all processes to finish before continuing
    wait(pre_process)

    print '''
####################################################################################
# Step 3: mark duplicate reads
####################################################################################
'''

    # mark duplicates
    mark_dups = [sample.mark_duplicates() for sample in samples]
    # wait for all processes to finish before continuing
    wait(mark_dups)

    print '''
####################################################################################
# Step 4: index bam files
####################################################################################
'''

    # index bam files
    index_bam = [sample.index_bam() for sample in samples]
    # wait for all processes to finish before continuing
    wait(index_bam)

    print '''
####################################################################################
# Step 5: call variants
####################################################################################
'''

    # size of chunks for variant call batching (based on hardware limitations)
    samples_per_chunk = 9

    # split sample list into chunks
    sample_chunks = [samples[i:i + samples_per_chunk] for i in xrange(0, len(samples), samples_per_chunk)]

    for chunk in sample_chunks:
        # call variants
        call_variants = [sample.call_variants() for sample in chunk]
        # wait for all processes to finish before continuing
        wait(call_variants)

    print '''
####################################################################################
# Step 6: create alignment summary file
####################################################################################
'''

    # TODO write aln summary

    print 'Pipeline complete!'



