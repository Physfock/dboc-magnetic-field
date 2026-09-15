Potential energy curves (PECs) for Diagonal Born–Oppenheimer Corrections (DBOC) as functions of magnetic field strength and internuclear distance for the H2, LiH, and HeH+ molecules. The repository contains a Podman container with modified ChronusQ code that writes the one- and two-electron matrix elements required for the subsequent correlation calculations in HDF5 format. This data is then reused in PySCF to calculate the DBOC within the approximation discussed in the manuscript: https://arxiv.org/pdf/2609.15327 

The following instructions are also available at: https://github.com/Physfock/cc_tesla 



Instructions:

Install Podman on Linux Refer to the official Podman installation guide: https://podman.io/docs/installation For Debian-based systems (e.g., Ubuntu): The podman package is available in Debian repositories: sudo apt-get -y install podman For Fedora-based systems: sudo dnf -y install podman
Download the container image Download the container image into your home directory (/home/<your_username>/) using the following link: https://app.filen.io/#/d/1000f5d4-52f0-433d-8fb6-3e28ec3a8dac%23uJ6zYPeFPCSkMtKWOBNBJAbJ4YTglNKO
Navigate to the download directory Open a terminal and navigate to your home directory (/home/<your_username>/), where you downloaded the .tar file in the previous step.
Load the image into Podman podman load -i /home/<your_username>/cc_tesla.tar
Verify the image was loaded podman images
Create a working directory Create a directory named Tesla in your home directory (/home/<your_username>/Tesla/), where your programs will be executed: mkdir Tesla 
Run the container podman run -it -v /home/<your_username>/Tesla/:/mnt/:Z cc_tesla bash This command starts the container and mounts your local /home/<your_username>/Tesla/ directory to /mnt/ inside the container. Any changes made to /mnt/ inside the container will immediately be reflected in your local Tesla directory, and vice versa. You are now inside the container's terminal. All necessary tools and software (e.g., nano, mc, Python, modified PySCF, and ChronusQ) are pre-installed.
Troubleshooting: SELinux permissions If you encounter permission issues or cannot see files in /home/<your_username>/Tesla/ from within the container, your system may be running SELinux in enforcing mode. Exit the container and temporarily disable SELinux with: sudo setenforce 0 
