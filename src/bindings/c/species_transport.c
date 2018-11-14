
#include "species_transport.h"

#include "string.h"

void initialise_transport_c(const char *MODEL, int *model_len);
void solve_transport_c(const char *MODEL, int *model_len);


void initialise_transport(const char *MODEL)
{
  int model_len = strlen(MODEL);

  initialise_transport_c(MODEL, &model_len);
}

void solve_transport(const char *MODEL)
{
  int model_len = strlen(MODEL);

  solve_transport_c(MODEL, &model_len);
}
