# IDEPA-XMBD

<div align=center><img src="./figs/IDEPA_log.svg" width="30%" height="30%" ></div>


We evaluated five state-of-the-art tools (RankComp v1/v2, PenDA, Peng, and Quantile) through classic computational (precision, Type one error control, parameter evaluation, robustness, and similarity) and functional (pathway enrichment and survival analysis) criteria. We also integrated these tools into a user-friendly tool kit, IDEPA-XMBD , to facilitate individualized DEAs in proteomics studies.

A pre-print describing the method is available at bioRxiv: [Application of personalized differential expression analysis in human cancer proteome](https://www.biorxiv.org/content/10.1101/2021.07.18.452812v2)

## Install
We use docker to encapsulate the command line version and the plotly version of IDEPA-XMBD separately.

### Cmd Version
Pull the docker image of IDEPA-XMBD:
```shell
docker pull ychlouie/idepa_cmd:latest
```

Create a docker container containing IDEPA-XMBD:
```shell
docker run -it ychlouie/idepa_cmd:latest
```

### Plotly Version

```shell

```

