# Instructions to Install and Run the Model

To install and run the model, follow these steps:

1. Open your terminal and run the following commands:
   ```
   sudo DOCKER_BUILDKIT=0 docker build -t feo-esmod-pypsa .
   docker run -it feo-esmod-pypsa:latest
   ```
2. In Visual Studio Code, add the Docker extension.

3. Go to the extensions and click on the "feo-esmod-pypsa" container.

4. Click on "Attach Visual Studio Code". This should open a new window.

5. Now, open the terminal and run the following command:
   ```
   snakemake c4 -j2 solve_all_networks
   ```
   This should run the tutorial for you.


That's it! You should now be able to run the model with ease.
