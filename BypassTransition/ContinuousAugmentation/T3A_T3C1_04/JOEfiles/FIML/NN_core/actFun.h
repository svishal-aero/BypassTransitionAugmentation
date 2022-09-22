double activate(double val, int actFun)
{
  switch(actFun)
  {
    case NN_LINEAR:  return val; break;
    case NN_SIGMOID: return 1.0 / (1.0 + exp(-val)); break;
    case NN_TANH:    return tanh(x); break;
    case NN_SWISH:   return val / (1.0 + exp(-val)); break;
  }
}

double d_activate(double val, int actFun)
{
  switch(actFun)
  {
    case NN_LINEAR:  return 1.0; break;
    case NN_SIGMOID: return exp(-val) / pow(1.0 + exp(-val), 2); break;
    case NN_TANH:    return 1.0 - tanh(x) * tanh(x); break;
    case NN_SWISH:   return 1.0 / (1.0 + exp(-val)) + val * exp(-val) / pow(1.0 + exp(-val), 2); break;
  }
}

adouble activate_AD(adouble val, int actFun)
{
  switch(actFun)
  {
    case NN_LINEAR:  return val; break;
    case NN_SIGMOID: return 1.0 / (1.0 + exp(-val)); break;
    case NN_TANH:    return tanh(x); break;
    case NN_SWISH:   return val / (1.0 + exp(-val)); break;
  }
}
