auto kmpi  = 139.57018;
auto kmp   = 938.2720813;
auto kmd   = 1875.612762;
auto kmt   = 2808.921112;
auto kmhe3 = 2808.39132;
auto kmal  = 3727.379378;
auto kme   = 0.5109989461;

/**
 * calculating dedx difference from reference to evaluated function value.
 * This function is for minimizing dedx function to calculate mass
 */
Double_t fddedx(Double_t *x, Double_t *p);

/**
 *  dE/dx function - Modified Bethe-Bloch formula
 */
Double_t fdedx(Double_t Z, Double_t m, Double_t *x, Double_t *p);

Double_t fdedx_pim (Double_t *x, Double_t *p) { return fdedx(-1, kmpi,  x, p); }
Double_t fdedx_pi  (Double_t *x, Double_t *p) { return fdedx( 1, kmpi,  x, p); }
Double_t fdedx_p   (Double_t *x, Double_t *p) { return fdedx( 1, kmp,   x, p); }
Double_t fdedx_d   (Double_t *x, Double_t *p) { return fdedx( 1, kmd,   x, p); }
Double_t fdedx_t   (Double_t *x, Double_t *p) { return fdedx( 1, kmt,   x, p); }
Double_t fdedx_he3 (Double_t *x, Double_t *p) { return fdedx( 2, kmhe3, x, p); }
Double_t fdedx_al  (Double_t *x, Double_t *p) { return fdedx( 2, kmal,  x, p); }

Double_t fddedx(Double_t *x, Double_t *p) {
  auto mass = x[0];
  // p[2] = mom
  // p[3] = Z
  // p[4] = dedx
  auto ZA = 17.2/37.6;
  auto I2 = 12.*p[3]+7.; I2 = I2*I2;
  auto pZ = p[2]*p[3];
  auto b2 = pZ*pZ/(pZ*pZ+mass*mass);
  auto g2 = 1./(1.-b2);
  auto Wx = 2*kme*b2*g2/((kme+mass)*(kme+mass)+2*kme*mass*(TMath::Sqrt(g2)-1));
  auto ddedx = p[4]-p[0]*ZA*p[3]*p[3]/b2*(0.5*TMath::Log(2*kme*b2*g2*Wx/I2)-b2-p[1]);
  return ddedx;
}

Double_t fdedx(Double_t Z, Double_t m, Double_t *x, Double_t *p) {
  auto ZA = 17.2/37.6;
  auto I2 = 12.*p[3]+7.; I2 = I2*I2;
  auto pZ = x[0]*Z;
  auto b2 = pZ*pZ/(pZ*pZ+m*m);
  auto g2 = 1./(1.-b2);
  auto Wx = 2*kme*b2*g2/((kme+m)*(kme+m)+2*kme*m*(TMath::Sqrt(g2)-1));
  auto dedx = p[0]*ZA*Z*Z/b2*(0.5*TMath::Log(2*kme*b2*g2*Wx/I2)-b2-p[1]);
  return dedx;
}
