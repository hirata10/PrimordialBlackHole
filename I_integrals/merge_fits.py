import sys
import os
import glob
from astropy.io import fits

# Get x from command-line arguments
x = int(sys.argv[1])

# If a results directory is provided as a second argument, use it.
if len(sys.argv) > 2:
    results_dir = sys.argv[2]
    file_pattern = os.path.join(results_dir, f'output_x{x}_omega*.fits')
else:
    file_pattern = f'output_x{x}_omega*.fits'

def get_omega_from_filename(filename):
    """Extracts the omega value from the filename 'output_x100_omega*.fits'."""
    # Assuming the filename is formatted as 'output_x100_omega<omega_value>.fits'
    return float(filename.split("omega")[1].split(".fits")[0])

def merge_fits_files(output_filename, file_pattern):
    # Get a list of all FITS files matching the pattern
    fits_files = glob.glob(file_pattern)
    
    if not fits_files:
        print("No FITS files found matching the pattern.")
        return

    # Sort the files based on the numerical value of omega
    fits_files = sorted(fits_files, key=get_omega_from_filename)
    
    print(f"Found {len(fits_files)} FITS files to merge.")
    
    # Initialize an HDUList to store the merged data
    hdulist = fits.HDUList()

    for i, filename in enumerate(fits_files):
        print(f"Merging file {i+1}: {filename}")       
        # Open the FITS file and copy the data and header before the file is closed
        with fits.open(filename) as hdul:
            data = hdul[0].data.copy()  # Copy the data to avoid closed file access
            header = hdul[0].header.copy()  # Copy the header

            # If it's the first file, use the Primary HDU
            if i == 0:
                hdulist.append(fits.PrimaryHDU(data=data, header=header))
            else:
                # For all subsequent files, append them as Image HDUs
                hdulist.append(fits.ImageHDU(data=data, header=header))

    # Write the merged FITS file
    try:
        hdulist.writeto(output_filename, overwrite=True)
        print(f"Merged file written to: {output_filename}")
    except Exception as e:
        print(f"Error writing merged file: {e}")

if __name__ == "__main__":
    # Output file name using the format with x
    output_filename = f"merged_output_x{x}.fits"   
    merge_fits_files(output_filename, file_pattern)
