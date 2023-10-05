<h1> RIME Sience Operation Center </h1>

<h2> Installation Guide:</h2>

- Install microk8s (following the official documentation at https://microk8s.io/) 
    ```
    sudo snap install microk8s --classic
    ```
- Install tekton pipelines (following the official documentation at https://tekton.dev/docs/installation/pipelines/) 
    ```
    microk8s kubectl apply --filename https://storage.googleapis.com/tekton-releases/pipeline/latest/release.yaml
    ```
- Install tekton triggers (following the official documentation at https://tekton.dev/docs/installation/triggers/)

    ```
    microk8s kubectl apply --filename \
    https://storage.googleapis.com/tekton-releases/triggers/latest/release.yaml

    microk8s kubectl apply --filename \
    https://storage.googleapis.com/tekton-releases/triggers/latest/interceptors.yaml
    ```

- Install tekton Dashboard (following official documentation at https://tekton.dev/docs/dashboard/install/)
    ```
    microk8s kubectl apply --filename https://storage.googleapis.com/tekton-releases/dashboard/latest/release.yaml
    ```

- Install ArgoCD (following the official documentation at https://argo-cd.readthedocs.io/en/stable/)
    ```
    microk8s kubectl create namespace argocd
    microk8s kubectl apply -n argocd -f https://raw.githubusercontent.com/argoproj/argo-cd/stable/manifests/install.yaml
    ```

- Access the ArgoCD dashboard and create two apps for the subfolders `Tekton` and `Deployments`, selecting `DIRECTORY RECURSE` as parameter and sync the two applications. This will deploy all the required configuration.

By accessing the tekton dashboard, all the tasks and pipelines should appear for the execution


<h2> Current status:</h2>

The current implementation uses a NFS server located on one node (to be configured). The file Deployments/Volumes/NFS_SC.yaml should be configured by changong `share` and `server` parameters with the personal environment.

The pipeline `telemetry2raw-raw2quicklook` takes the telemetry passed as input, that should be located in `data` inside the mounted storage. The pipeline will provide the raw file, the quicklook products and the logs of the steps in the dedicated folders that will be created if not present.

<h2> Future implementations:</h2>

- Better management of the volumes: create a dedicated PVC for every pipeline
- Better resource management: the tm2raw.py script is using all the available computational resources. It may be better to limit the utilization of the script in order to ensure a good parallelization since the script should process only one telemetry at execution
- Create a trigger