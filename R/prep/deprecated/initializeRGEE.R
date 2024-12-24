### try google earth engine R for the 100th time

install.packages(c("remotes", "googledrive"))
remotes::install_github("r-spatial/rgee")

library(rgee)

# Get the username
HOME <- Sys.getenv("HOME")

# 1. Install miniconda
reticulate::install_miniconda()

# 2. Install Google Cloud SDK
system("curl -sSL https://sdk.cloud.google.com | bash")

# 3 Set global parameters
Sys.setenv("RETICULATE_PYTHON" = sprintf("%s/.local/share/r-miniconda/bin/python3", HOME))
Sys.setenv("EARTHENGINE_GCLOUD" = sprintf("%s/google-cloud-sdk/bin/", HOME))


ee$Authenticate(auth_mode='notebook')

# Initialize - this will connect to a project. You should always call this
# before working with rgee. It is IMPORTANT THAT YOU SPECIFY A PROJECT using
# the project parameter. If you forget what project IDs you have access to, find them
# here: console.cloud.google.com/project
ee$Initialize(project='ee-jonastrepel')  # <-- EDIT THIS FOR YOUR PROJECT

# Optionally make a request to verify you are connected.
ee$String('Hello from the Earth Engine servers!')$getInfo()

### It worked! 
