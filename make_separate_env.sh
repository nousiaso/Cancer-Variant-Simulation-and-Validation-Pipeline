#!/bin/bash

# First, ensure channels are in correct order
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Remove existing environments if they exist
conda deactivate
for env in bam_sim_base bam_sim_samtools bam_sim_somatosim bam_sim_snake; do
    conda env remove -n $env -y
done

# Create base environment with Python packages
echo "Creating base environment..."
mamba create -n bam_sim_base -y python=3.9 \
    pip \
    numpy \
    pandas \
    scipy \
    scikit-learn \
    matplotlib=3.9.2 \
    seaborn \
    pytest \
    pysam=0.19.1 \
    biopython \
    requests \
    click \
    tqdm

# Create SomatoSim environment with specific versions
echo "Creating SomatoSim environment..."
mamba create -n bam_sim_somatosim -y \
    python=3.7 \
    cycler=0.10.0 \
    cython=0.29.15 \
    kiwisolver=1.1.0 \
    matplotlib=3.2.0 \
    numpy=1.18.5 \
    pandas=1.0.2 \
    pyparsing=2.4.6 \
    pysam=0.15.4 \
    python-dateutil=2.8.1 \
    pytz=2019.3 \
    six=1.14.0 \
    samtools \
    scikit-learn=1.0.2 \
    scipy=1.7.3 \
    seaborn=0.12.2

# Activate SomatoSim environment and install SomatoSim
conda activate bam_sim_somatosim
cd SomatoSim
python -m pip install .
cd ..
conda deactivate

# Create Samtools environment with specific version
echo "Creating Samtools environment..."
mamba create -n bam_sim_samtools -y \
    samtools=1.9 \
    htslib=1.9

# Create Snakemake environment
echo "Creating Snakemake environment..."
mamba create -n bam_sim_snake -y \
    python=3.9 \
    snakemake-minimal=7.32.4

# Create wrapper scripts
mkdir -p ~/bin

# Add SomatoSim wrapper script
cat > ~/bin/somatosim_wrapper.sh << 'EOF'
#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bam_sim_somatosim
somatosim "$@"
conda deactivate
EOF

# Add Samtools wrapper script
cat > ~/bin/samtools_wrapper.sh << 'EOF'
#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bam_sim_samtools
samtools "$@"
conda deactivate
EOF

# Add Snakemake wrapper script
cat > ~/bin/snakemake_wrapper.sh << 'EOF'
#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bam_sim_snake
snakemake "$@"
conda deactivate
EOF

# Make wrapper scripts executable
chmod +x ~/bin/somatosim_wrapper.sh
chmod +x ~/bin/samtools_wrapper.sh
chmod +x ~/bin/snakemake_wrapper.sh

echo "Installation complete. Created environments:"
echo "- bam_sim_base (Python packages)"
echo "- bam_sim_somatosim (SomatoSim)"
echo "- bam_sim_samtools (Samtools)"
echo "- bam_sim_snake (Snakemake workflow manager)"

echo -e "\nWrapper scripts created in ~/bin:"
echo "- somatosim_wrapper.sh"
echo "- samtools_wrapper.sh"
echo "- snakemake_wrapper.sh"

echo -e "\nAdd this line to your ~/.bashrc or ~/.zshrc:"
echo 'export PATH="$HOME/bin:$PATH"'

echo -e "\nUsage examples:"
echo "somatosim_wrapper.sh --help"
echo "samtools_wrapper.sh --version"
echo "snakemake_wrapper.sh --version"
