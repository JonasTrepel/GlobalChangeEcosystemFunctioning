### try google earth engine R for the 100th time

# install.packages(c("remotes", "googledrive"))
# remotes::install_github("r-spatial/rgee")

## macbook path /opt/anaconda3/envs/rgee_env/bin
## GIS04 server: "C:\\Users\\au713983\\.conda\\envs\\rgee_env" 
## AU laptop: "C:\\Users\\au713983\\AppData\\Local\\anaconda3\\envs\\rgee_env"
rgee_env_dir <- c("/opt/anaconda3/envs/rgee_env/bin")

library(rgee)

reticulate::use_python(rgee_env_dir, required=T)
rgee::ee_install_set_pyenv(
  py_path = rgee_env_dir,
  py_env = 'rgee_env'
)

#Set system environment variables
Sys.setenv(RETICULATE_PYTHON = rgee_env_dir)
Sys.setenv(EARTHENGINE_PYTHON = rgee_env_dir)

#### AUTHENTICATE AGAIN#####
ee_clean_user_credentials()
ee$Authenticate(auth_mode='notebook')

# Initialize - this will connect to a project. You should always call this
# before working with rgee. It is IMPORTANT THAT YOU SPECIFY A PROJECT using
# the project parameter. If you forget what project IDs you have access to, find them
# here: console.cloud.google.com/project
ee$Initialize(project='ee-jonastrepel')  # <-- EDIT THIS FOR YOUR PROJECT

# Optionally make a request to verify you are connected.
ee$String('Hello from the Earth Engine servers!')$getInfo()

# Create visualization parameters
viz_params <- list(min = 0, max = 1000)

# Add the layer to the map
Map$centerObject(elev, 6)  # Centers the map on the elevation data
Map$addLayer(elev, list(min = 0, max = 7000))



