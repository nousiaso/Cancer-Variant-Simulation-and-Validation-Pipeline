#!/bin/bash

# First, ensure channels are in correct order
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Remove existing environments if they exist
conda deactivate
for env in bam_sim_base bam_sim_art bam_sim_samtools bam_sim_bwa bam_sim_snake bam_sim_somatosim; do
    conda env remove -n $env -y
done

# Create base environment with Python packages
echo "Creating base environment..."
mamba create -n bam_sim_base -y python=3.9 \
    pip \
    numpy \
    pandas \
    scipyq \
    requests \
    biopython \
    pytest \
    click \
    tqdm \
    matplotlib \
    seaborn \
    pysam=0.19.1

# Create SomatoSim environment with specific versions
echo "Creating SomatoSim environment..."
mamba create -n bam_sim_somatosim -y \
    python=3 \
    cycler=0.10.0 \
    Cython=0.29.15 \
    kiwisolver=1.1.0 \
    matplotlib=3.2.0 \
    numpy=1.18 \
    pandas=1.0.2 \
    pyparsing=2.4.6 \
    pysam=0.15.4 \
    python-dateutil=2.8.1 \
    pytz=2019.3 \
    six=1.14.0 \
    samtools

# Activate SomatoSim environment and install SomatoSim
conda activate bam_sim_somatosim
cd SomatoSim
python -m pip install .
cd ..
conda deactivate

# [Previous environment creations remain the same...]
# Create ART environment
echo "Creating ART environment..."
mamba create -n bam_sim_art -y art

# Create Samtools environment
echo "Creating Samtools environment..."
mamba create -n bam_sim_samtools -y \
    libdeflate=1.17 \
    samtools=1.16

# Create BWA environment
echo "Creating BWA environment..."
mamba create -n bam_sim_bwa -y bwa=0.7.17

# Create Snakemake environment
echo "Creating Snakemake environment..."
mamba create -n bam_sim_snake -y \
    python=3.9 \
    snakemake-minimal=7.32.4

# [Previous test functions remain the same...]

# Add SomatoSim wrapper script
cat > ~/bin/somatosim_wrapper.sh << 'EOF'
#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bam_sim_somatosim
somatosim "$@"
conda deactivate
EOF

# Make SomatoSim wrapper script executable
chmod +x ~/bin/somatosim_wrapper.sh

# [Previous wrapper script creations remain the same...]
cat > ~/bin/art_wrapper.sh << 'EOF'
#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bam_sim_art
art_illumina "$@"
conda deactivate
EOF

cat > ~/bin/samtools_wrapper.sh << 'EOF'
#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bam_sim_samtools
samtools "$@"
conda deactivate
EOF

cat > ~/bin/bwa_wrapper.sh << 'EOF'
#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bam_sim_bwa
bwa "$@"
conda deactivate
EOF

cat > ~/bin/snakemake_wrapper.sh << 'EOF'
#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bam_sim_snake
snakemake "$@"
conda deactivate
EOF

# Make all wrapper scripts executable
chmod +x ~/bin/art_wrapper.sh
chmod +x ~/bin/samtools_wrapper.sh
chmod +x ~/bin/bwa_wrapper.sh
chmod +x ~/bin/snakemake_wrapper.sh

echo "Installation complete. Created environments:"
echo "- bam_sim_base (Python packages)"
echo "- bam_sim_somatosim (SomatoSim)"
echo "- bam_sim_art (ART)"
echo "- bam_sim_samtools (Samtools)"
echo "- bam_sim_bwa (BWA)"
echo "- bam_sim_snake (Snakemake)"

echo -e "\nWrapper scripts created in ~/bin:"
echo "- somatosim_wrapper.sh"
echo "- art_wrapper.sh"
echo "- samtools_wrapper.sh"
echo "- bwa_wrapper.sh"
echo "- snakemake_wrapper.sh"

echo -e "\nAdd this line to your ~/.bashrc or ~/.zshrc:"
echo 'export PATH="$HOME/bin:$PATH"'

echo -e "\nUsage examples:"
echo "somatosim_wrapper.sh --help"
echo "art_wrapper.sh --version"
echo "samtools_wrapper.sh --version"
echo "bwa_wrapper.sh"
echo "snakemake_wrapper.sh --version"
