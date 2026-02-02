#!/bin/bash
# =========================================================================================
# SCRIPT: build_cellranger.sh
# =========================================================================================
# Purpose: Build Cell Ranger Singularity/Apptainer image for HPC or local environments
#
# What it does:
#   1. Downloads Cell Ranger tarball from 10x Genomics (if not present)
#   2. Builds Singularity image using cellranger.def definition file
#   3. Validates image integrity by checking Cell Ranger version
#   4. Cleans up large tarball after successful build
#   5. Logs all output to timestamped log file
#
# Prerequisites:
#   - Singularity/Apptainer installed (via module system or system-wide)
#   - Sufficient disk space (~20GB for tarball + ~10GB for final image)
#   - cellranger.def file in same directory as this script
#   - Internet connection for tarball download (if not cached)
#
# Output:
#   - cellranger_10.0.0.sif: Ready-to-use Singularity image
#   - build_log_YYYYMMDD_HHMMSS.txt: Complete build log
#
# Usage on HPC (with root/fakeroot access):
#   bash build_cellranger.sh
#
# Typical runtime: 10-20 minutes (depends on download speed and I/O)
# =========================================================================================

# =========================================================================================
# WINDOWS/WSL BUILD INSTRUCTIONS (No HPC Root Access)
# =========================================================================================
# If you DON'T have root access on HPC, build the image on Windows using WSL:
#
# Step 1: Install WSL and Ubuntu
#   Open Windows PowerShell as Administrator:
#     wsl --install -d Ubuntu
#   Reboot computer when prompted
#
# Step 2: Install Apptainer in WSL
#   Open WSL Ubuntu terminal:
#     sudo apt update && sudo apt upgrade -y
#     sudo apt install -y software-properties-common
#     sudo add-apt-repository -y ppa:apptainer/ppa
#     sudo apt update
#     sudo apt install -y apptainer
#
# Step 3: Navigate to build directory
#   Access your Windows files from WSL:
#     cd "/mnt/c/Users/YOUR_USERNAME/path/to/build/directory"
#   Example:
#     cd "/mnt/c/Users/kailasamms/Desktop/build"
#
# Step 4: Verify files are present
#     ls -lh
#   Should see: cellranger.def, cellranger-10.0.0.tar.gz (or will download)
#
# Step 5: Build the image
#     sudo apptainer build cellranger_10.0.0.sif cellranger.def
#   Note: This may take 15-30 minutes
#
# Step 6: Verify image
#     singularity exec cellranger_10.0.0.sif cellranger --version
#   Should output: cellranger cellranger-10.0.0
#
# Step 7: Upload to HPC
#   - Copy cellranger_10.0.0.sif to your HPC singularity_cache directory
#   - Use scp, rsync, or your institution's data transfer tool
#   Example:
#     scp cellranger_10.0.0.sif username@hpc.edu:/path/to/singularity_cache/
#
# Troubleshooting WSL:
#   - "Permission denied" → Use sudo for apptainer commands
#   - "No space left" → WSL default disk is small, increase in .wslconfig
#   - "Can't access Windows files" → Use /mnt/c/ path, not C:\
#   - Build fails → Check cellranger.def file is present and valid
# =========================================================================================

set -e  # Exit immediately if any command fails
set -u  # Exit if undefined variable is used
set -o pipefail  # Exit if any command in a pipeline fails

# =========================================================================================
# CONFIGURATION VARIABLES
# =========================================================================================

# Cell Ranger version and download details
VERSION="10.0.0"
TAR_FILE="cellranger-${VERSION}.tar.gz"
IMG_NAME="cellranger_${VERSION}.sif"
IMG_SCRIPT="cellranger.def"

# Download URL from 10x Genomics (with signed authentication)
# Note: This URL may expire - get fresh URL from: https://www.10xgenomics.com/support/software/cell-ranger/downloads
DOWNLOAD_URL="https://cf.10xgenomics.com/releases/cell-exp/cellranger-10.0.0.tar.gz?Expires=1769765821&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=YerBr4hrgtShcJX-IejhWom7wwOStFD-GjakodZ7UlrMLAS9mfQeLmfEx3TNbQG2op7DSJnZW3qbUUxQdMu7apMyeQz5iK36ZsJcrk5TAeMZX8y2k9oxZyIB6xowDjzMMIICd14oOZJcQNpsSuClbvrPyThGQW3JHa4MBrJJKp~-r1Kuq07nOK9ggSg6R9sNd04RLBlpBAjd2Y6~kxcGnbsuD8jIAmBPGpGSkVZbhsjIIfRWdTSEs-nObwfU6w-0BXwsTq~0~saVfVTLIfB03lg3UXu7sE3FHGJ4gn7rN9lfrjqsl6ezBuIJix3lPa34lUX8CAC2GSCSTBj8S3H2Lg__"

# Logging
LOG_FILE="build_log_$(date +%Y%m%d_%H%M%S).txt"

# =========================================================================================
# INITIALIZATION
# =========================================================================================

# Start logging everything to both screen and log file
exec &> >(tee -a "${LOG_FILE}")

echo "================================================================================"
echo "Cell Ranger ${VERSION} Singularity Image Build"
echo "================================================================================"
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo "Log file: ${LOG_FILE}"
echo "================================================================================"
echo ""

# Move to script directory (handles being called from other locations)
cd "$(dirname "$(readlink -f "$0")")" || {
    echo "ERROR: Cannot change to script directory"
    exit 1
}

# =========================================================================================
# ENVIRONMENT SETUP
# =========================================================================================

echo "Step 1: Setting up environment..."
echo ""

# Try to load Singularity/Apptainer module (HPC environments)
# The || true prevents script from exiting if module command fails
if command -v module &> /dev/null; then
    echo "Module system detected. Attempting to load Singularity/Apptainer..."
    # Try multiple common module names
    module load singularity-apptainer/1.1.8 2>/dev/null || \
    module load singularity 2>/dev/null || \
    module load apptainer 2>/dev/null || \
    echo "  → No module loaded (assuming system-wide installation)"
else
    echo "No module system detected (assuming system-wide Singularity/Apptainer)"
fi

# Verify Singularity/Apptainer is available
if ! command -v singularity &> /dev/null && ! command -v apptainer &> /dev/null; then
    echo ""
    echo "ERROR: Neither 'singularity' nor 'apptainer' command found!"
    echo ""
    echo "Please install Singularity/Apptainer or load the appropriate module."
    echo "For installation instructions, see:"
    echo "  - Singularity: https://docs.sylabs.io/guides/latest/user-guide/"
    echo "  - Apptainer: https://apptainer.org/docs/user/latest/"
    exit 1
fi

# Display version
echo ""
if command -v singularity &> /dev/null; then
    echo "✓ Singularity version: $(singularity --version)"
elif command -v apptainer &> /dev/null; then
    echo "✓ Apptainer version: $(apptainer --version)"
fi

# Set temporary directory to current location (avoids filling up /tmp)
export SINGULARITY_TMPDIR=$PWD
export APPTAINER_TMPDIR=$PWD
echo "✓ Temporary directory: ${SINGULARITY_TMPDIR}"
echo ""

# =========================================================================================
# PREREQUISITE CHECKS
# =========================================================================================

echo "Step 2: Checking prerequisites..."
echo ""

# Check for definition file
if [ ! -f "${IMG_SCRIPT}" ]; then
    echo "ERROR: Definition file not found: ${IMG_SCRIPT}"
    echo ""
    echo "Please ensure ${IMG_SCRIPT} is in the same directory as this script."
    echo "Current directory: $(pwd)"
    echo "Files present:"
    ls -lh
    exit 1
fi
echo "✓ Definition file found: ${IMG_SCRIPT}"

# Check available disk space
AVAILABLE_SPACE=$(df -BG . | awk 'NR==2 {print $4}' | sed 's/G//')
REQUIRED_SPACE=30  # GB needed for safe build
if [ "${AVAILABLE_SPACE}" -lt "${REQUIRED_SPACE}" ]; then
    echo ""
    echo "WARNING: Low disk space detected!"
    echo "  Available: ${AVAILABLE_SPACE}GB"
    echo "  Recommended: ${REQUIRED_SPACE}GB"
    echo ""
    read -p "Continue anyway? (y/n): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Build cancelled by user."
        exit 1
    fi
else
    echo "✓ Sufficient disk space: ${AVAILABLE_SPACE}GB available"
fi
echo ""

# =========================================================================================
# TARBALL DOWNLOAD
# =========================================================================================

echo "Step 3: Obtaining Cell Ranger tarball..."
echo ""

if [ -f "${TAR_FILE}" ]; then
    echo "✓ Tarball already exists: ${TAR_FILE}"
    TARBALL_SIZE=$(du -h "${TAR_FILE}" | cut -f1)
    echo "  Size: ${TARBALL_SIZE}"
    echo "  Skipping download."
else
    echo "Tarball not found. Downloading from 10x Genomics..."
    echo "Source: ${DOWNLOAD_URL:0:80}..."
    echo ""
    
    # Check if wget is available
    if ! command -v wget &> /dev/null; then
        echo "ERROR: wget not found. Please install wget or manually download:"
        echo "${DOWNLOAD_URL}"
        echo "Save as: ${TAR_FILE}"
        exit 1
    fi
    
    # Download with progress bar and resume capability
    # -c: Continue partial downloads
    # -O: Output filename
    # --progress=bar: Show progress
    wget -c -O "${TAR_FILE}" "${DOWNLOAD_URL}" || {
        echo ""
        echo "ERROR: Download failed!"
        echo ""
        echo "Possible causes:"
        echo "  1. Download URL expired (10x uses signed URLs with expiration)"
        echo "  2. Network connectivity issues"
        echo "  3. Firewall blocking downloads"
        echo ""
        echo "Solution: Get fresh download URL from:"
        echo "  https://www.10xgenomics.com/support/software/cell-ranger/downloads"
        echo ""
        echo "Update DOWNLOAD_URL variable in this script with the new URL."
        exit 1
    }
    
    echo ""
    echo "✓ Download completed: ${TAR_FILE}"
fi
echo ""

# =========================================================================================
# IMAGE BUILD
# =========================================================================================

echo "Step 4: Building Singularity image..."
echo ""
echo "This may take 10-20 minutes depending on system I/O..."
echo "Image: ${IMG_NAME}"
echo "Definition: ${IMG_SCRIPT}"
echo ""

# Record build start time
BUILD_START=$(date +%s)

# Build command explanation:
# --fakeroot: Build without root (uses user namespaces on HPC)
# ${IMG_NAME}: Output image filename
# ${IMG_SCRIPT}: Input definition file
#
# Note: If --fakeroot fails (not supported), try:
#   sudo singularity build ${IMG_NAME} ${IMG_SCRIPT}

singularity build --fakeroot "${IMG_NAME}" "${IMG_SCRIPT}" || {
    BUILD_EXIT_CODE=$?
    echo ""
    echo "================================================================================"
    echo "BUILD FAILED"
    echo "================================================================================"
    echo ""
    echo "Exit code: ${BUILD_EXIT_CODE}"
    echo ""
    echo "Common causes:"
    echo "  1. --fakeroot not supported → Try with sudo (if available)"
    echo "  2. Insufficient space → Free up disk space"
    echo "  3. Definition file errors → Check ${IMG_SCRIPT} syntax"
    echo "  4. Network issues → Check if tarball extraction worked"
    echo ""
    echo "Check full log: ${LOG_FILE}"
    echo "Keeping ${TAR_FILE} for debugging."
    echo ""
    echo "================================================================================"
    exit 1
}

# Calculate build time
BUILD_END=$(date +%s)
BUILD_TIME=$((BUILD_END - BUILD_START))
BUILD_MINUTES=$((BUILD_TIME / 60))
BUILD_SECONDS=$((BUILD_TIME % 60))

echo ""
echo "✓ Build completed successfully"
echo "  Time: ${BUILD_MINUTES}m ${BUILD_SECONDS}s"
echo "  Image: ${IMG_NAME}"
IMG_SIZE=$(du -h "${IMG_NAME}" | cut -f1)
echo "  Size: ${IMG_SIZE}"
echo ""

# =========================================================================================
# IMAGE VALIDATION
# =========================================================================================

echo "Step 5: Validating image integrity..."
echo ""
echo "Running: singularity exec ${IMG_NAME} cellranger --version"
echo ""

# Critical validation: Verify Cell Ranger works in the image
if singularity exec "${IMG_NAME}" cellranger --version; then
    echo ""
    echo "================================================================================"
    echo "✓ VALIDATION PASSED: Cell Ranger is functional"
    echo "================================================================================"
    echo ""
else
    echo ""
    echo "================================================================================"
    echo "✗ VALIDATION FAILED: Cell Ranger version check failed"
    echo "================================================================================"
    echo ""
    echo "The image built successfully but Cell Ranger cannot execute."
    echo "This usually indicates:"
    echo "  1. Missing dependencies in definition file"
    echo "  2. Incorrect PATH settings"
    echo "  3. Library incompatibilities"
    echo ""
    echo "Action: Keeping ${TAR_FILE} for debugging."
    echo "Check build log: ${LOG_FILE}"
    echo ""
    echo "DO NOT use this image for production work!"
    echo "================================================================================"
    exit 1
fi

# =========================================================================================
# CLEANUP
# =========================================================================================

echo "Step 6: Cleaning up..."
echo ""

# Only remove tarball if validation passed
echo "Removing large tarball: ${TAR_FILE}"
rm -f "${TAR_FILE}"
echo "✓ Cleanup completed"
echo ""

# =========================================================================================
# COMPLETION SUMMARY
# =========================================================================================

echo "================================================================================"
echo "BUILD COMPLETED SUCCESSFULLY"
echo "================================================================================"
echo ""
echo "Image details:"
echo "  Filename: ${IMG_NAME}"
echo "  Size: ${IMG_SIZE}"
echo "  Version: Cell Ranger ${VERSION}"
echo ""
echo "Next steps:"
echo "  1. Test the image:"
echo "     singularity exec ${IMG_NAME} cellranger --help"
echo ""
echo "  2. Move to singularity cache directory (if on HPC):"
echo "     mv ${IMG_NAME} /path/to/singularity_cache/"
echo ""
echo "  3. Use in Nextflow pipeline:"
echo "     Update singularity_cache path in nextflow.config"
echo ""
echo "Build log saved: ${LOG_FILE}"
echo "Finished: $(date)"
echo "================================================================================"

exit 0

# =========================================================================================
# TROUBLESHOOTING REFERENCE
# =========================================================================================
#
# Problem: "fakeroot: command not found" or build fails with --fakeroot
# Solution: System doesn't support fakeroot. Try:
#           sudo singularity build cellranger_10.0.0.sif cellranger.def
#           Or build on Windows/WSL (see instructions at top of file)
#
# Problem: "No space left on device"
# Solution: Free up space or set SINGULARITY_TMPDIR to larger partition:
#           export SINGULARITY_TMPDIR=/scratch/$USER
#           bash build_cellranger.sh
#
# Problem: Download URL expired (403 Forbidden)
# Solution: Get fresh URL from 10x Genomics website:
#           1. Visit: https://www.10xgenomics.com/support/software/cell-ranger/downloads
#           2. Right-click download link → Copy link address
#           3. Update DOWNLOAD_URL variable in this script
#
# Problem: "cellranger: command not found" after build
# Solution: Definition file issue. Verify:
#           1. Tarball extracted correctly
#           2. PATH set correctly in definition file
#           3. Cell Ranger binaries are executable
#
# Problem: Validation passes but pipeline fails
# Solution: Test Cell Ranger commands manually:
#           singularity exec cellranger_10.0.0.sif cellranger testrun --id=tiny
#           Check for missing dependencies or library issues
#
# =========================================================================================
