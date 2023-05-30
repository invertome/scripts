import subprocess
import sys
import os

def install(package, local=False):
    if local:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "--user", package])
    else:
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])

def is_module_installed(module_name):
    try:
        __import__(module_name)
        return True
    except ImportError:
        return False

def check_and_install(module_name, package_name=None):
    if package_name is None:
        package_name = module_name

    if not is_module_installed(module_name):
        print(f"{package_name} is not installed.")
        response = input("Do you want to install it? (Y/n): ").strip()
        if response.lower() in {"", "y", "yes"}:
            try:
                install(package_name)
                print(f"{package_name} has been installed successfully.")
            except subprocess.CalledProcessError:
                print("Installation failed. Trying to install locally.")
                try:
                    install(package_name, local=True)
                    print(f"{package_name} has been installed locally successfully.")
                except subprocess.CalledProcessError:
                    print(f"Failed to install {package_name} locally.")

def install_pip():
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pip"])
        print("pip has been installed successfully.")
    except subprocess.CalledProcessError:
        print("Installation failed. Trying to install locally.")
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", "--user", "pip"])
            print("pip has been installed locally successfully.")
        except subprocess.CalledProcessError:
            print("Failed to install pip locally.")

# First check for pip
if not is_module_installed('pip'):
    print("pip is not installed.")
    response = input("Do you want to install it? (Y/n): ").strip()
    if response.lower() in {"", "y", "yes"}:
        install_pip()
    else:
        print("pip is needed for this script to check and install other dependencies.")
        sys.exit(1)

dependencies = {
    "numpy": "numpy",
    "pandas": "pandas",
    "matplotlib": "matplotlib",
    "Bio": "biopython",
    "logomaker": "logomaker",
    "gecos": "gecos"
}

for module_name, package_name in dependencies.items():
    check_and_install(module_name, package_name)
