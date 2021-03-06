# SparseSensors.jl
   
   This repository is an implementation of the core sparse sensor placement with QR factorization and cost-constrained QR factorization algorithms from Manohar, _et al.,_ "Data-Driven Sparse Sensor Placcement for Reconstruction", and other papers, in Julia. This is a hobbyist port of the fantastic Python library [pysensors](https://github.com/dynamicslab/pysensors) in Julia. 
   
   All collaborations and contributions are welcome.
   
   
   
   ## Installation
   To install, use Pkg. From the Julia REPL, press ] to enter Pkg-mode and run
   ```julia
   pkg> add SparseSensors
   ```
   
   ## Example
   
   ```julia
   using SparseSensors
   using LinearAlgebra
   using Gadfly
   using DataFrames
   
   #Setup the experiment
   r = 11; # Number of basis modes
   n = 1000;
   x = collect(0.0:1/n:1.0);
   vde_basis = VandermondeBasis(x,r);
   Ψ = vde_basis.Ψ; #Get the vandermonde basis matrix from the Basis struct
   
   #Make design matrix
   X = copy(transpose(Ψ));
   n_samples,n_features = size(X);
   
   #Setup QR pivot sensor placement algorithm
   qr_pivot = QRPivot(X);
   fit(qr_pivot);
   pivots = qr_pivot.pivots;
   
   #Select the top 15 sensor locations
   f = abs.(x.^2 .- 0.5);
   selected_sensors = get_sensors(pivots,15);
   x_sensed = x[selected_sensors];
   y_sensed = f[selected_sensors];
   
   #Ground truth
   df_true = DataFrame();
   df_true[!,"x_true"] = x
   df_true[!,"y_true"] = f
   
   #Sensed
   df = DataFrame()
   df[!,"x_sensed"] = x_sensed;
   df[!,"y_sensed"] = y_sensed;
   
   #Plot the results
   p1 = plot(df,
       layer(x=:x_sensed,y=:y_sensed,color=["Optimized Sensors"]),
       layer(df_true,x=:x_true,y=:y_true,Geom.line,Geom.point,color=["True Function"]),
       Guide.xlabel("x"),Guide.ylabel("y"))
   ```
   
   ![](example.svg)
   
   
   ## Dependencies
   
   julia >v1.6.
   LinearAlgebra, Gadfly and DataFrames for plotting.
   
   ## Todo:
   
   - Implement high level SSPOR and SSPOC interfaces
   - Implement high level basis representation wrapper for SVD, etc.
   - Set up compatibility with JuliaML/MLJ.jl/ScikitLearn.jl
   
   References
   ------------
   
   -  Manohar, Krithika, Bingni W. Brunton, J. Nathan Kutz, and Steven L. Brunton.
      "Data-driven sparse sensor placement for reconstruction: Demonstrating the
      benefits of exploiting known patterns."
      IEEE Control Systems Magazine 38, no. 3 (2018): 63-86.
      `[DOI] <https://doi.org/10.1109/MCS.2018.2810460>`
   
   -  Clark, Emily, Travis Askham, Steven L. Brunton, and J. Nathan Kutz.
      "Greedy sensor placement with cost constraints." IEEE Sensors Journal 19, no. 7
      (2018): 2642-2656.
      `[DOI] <https://doi.org/10.1109/JSEN.2018.2887044>`
   
   -  de Silva, Brian M., Krithika Manohar, Emily Clark, Bingni W. Brunton,
      Steven L. Brunton, J. Nathan Kutz.
      "PySensors: A Python package for sparse sensor placement."
      arXiv preprint arXiv:2102.13476 (2021). `[arXiv] <https://arxiv.org/abs/2102.13476>`
