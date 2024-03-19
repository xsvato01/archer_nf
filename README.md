# Nextflow Archer pipeline for CMBG @2024

**Make sure to get the correct config parameters**. This is a private repository, you will need to create secrets to pull this repository by Nextflow:

append to `~/.nextflow/scm`:
`providers {
    github {
        user = 'UCO'
        password = 'YOUPERSONALTOKEN'
        }
    }`

## If your docker image is private and hosted on Gitlab Muni:

Run docker login registry.gitlab.ics.muni.cz:443 -u <username> -p <token> on an ubuntu machine. Find the ~/.docker/config.json file. Run base64 config.json and copy the output into dockerconfig-secret.yaml (in git repo). Finally run kubectl --namespace [your new namespace] apply -f dockerconfig-secret.yaml to add the secret to the kuba cluster.
This secret has to be added to nextflow.config:
`process {
    pod = [[imagePullSecret:'gitlab-svaton-secret'], ...]
    }`

## Samplesheet

There is a online CHANGE THIS [samplesheet](https://docs.google.com/spreadsheets/d/1WOktQDMH13d_zr0g8be1BqysMd7TIyy4LrzPywfvjhA/edit#gid=0) to fill.
Make sure to run `getSamplesheet.sh` to download and convert the samplesheet to input Nextlow parameters.

@Jan Svaton 2024
