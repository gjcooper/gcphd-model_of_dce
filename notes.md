# Notes from meeting Scott and guy
We use Multivariate normals, which are then transformed to positive, scaling

## Possible options
Stick-breaking - bad priors
Dirichlet - exactly how to do this though

## Plan
For each of the 7 models have a count parameter
enter ll func
Instead of calc log like for all 7 models, sample from dirichlet based on count parameters (exp'd)
This gives the winning model for this particle. Then generate the log likelihood for this winning model.

Start points for counts -> [1,1,1,1,1,1,1]
prior for each count -> N(1, 2)

ll() {
  x = exp(x)
  model = rdirichlet(1, counts)
  like = ddensity(model)
  }

  # Can we store the model that wins?
