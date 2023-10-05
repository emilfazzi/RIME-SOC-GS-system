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

- Access on ArgoCD 