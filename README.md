# Instructions to Install and Run the Model

To install and run the model, follow these steps:

1. Open your terminal and run the following commands:

   sudo DOCKER_BUILDKIT=0 docker build -t pypsa-earth .
   docker run -it pypsa-earth:latest
   ```

2. In Visual Studio Code, add the Docker extension.

3. Go to the extensions and click on the "pypsa-earth" container.

4. Click on "Attach Visual Studio Code". This should open a new window.

   ```
   snakemake c4 -j2 solve_all_networks
   ```

   This should run the tutorial for you.


That's it! You should now be able to install and run the model with ease.
