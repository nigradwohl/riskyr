## Interfacing table calculation with the riskyr function (development):



riskyr <- function(scen.lbl = "",  ## WAS: txt$scen.lbl,
                   scen.lng = txt$scen.lng,
                   scen.txt = txt$scen.txt, popu.lbl = txt$popu.lbl,
                   cond.lbl = txt$cond.lbl,
                   cond.true.lbl = txt$cond.true.lbl, cond.false.lbl = txt$cond.false.lbl,
                   dec.lbl = txt$dec.lbl,
                   dec.pos.lbl = txt$dec.pos.lbl, dec.neg.lbl = txt$dec.neg.lbl,
                   hi.lbl = txt$hi.lbl, mi.lbl = txt$mi.lbl,
                   fa.lbl = txt$fa.lbl, cr.lbl = txt$cr.lbl,
                   prev = NA,
                   sens = NA,
                   spec = NA, fart = NA,
                   N = NA,  ## WAS: freq$N,
                   hi = NA, mi = NA,
                   fa = NA, cr = NA,
                   scen.src = txt$scen.src,
                   scen.apa = txt$scen.apa) {

  ## Note:  this interface is specific to the diagnostic case (alias: riskyr.diagnostic):


  ## Matrix for condition-decision case:
  cnd_dec <- matrix()


  ## (0): Initialize some stuff: ------
  num_mat <- matrix()

  object <- list(scen.lbl = scen.lbl, scen.lng = scen.lng, scen.txt = scen.txt,
                 popu.lbl = popu.lbl, cond.lbl = cond.lbl,
                 cond.true.lbl = cond.true.lbl, cond.false.lbl = cond.false.lbl,
                 dec.lbl = dec.lbl, dec.pos.lbl = dec.pos.lbl, dec.neg.lbl = dec.neg.lbl,
                 hi.lbl = hi.lbl, mi.lbl = mi.lbl, fa.lbl = fa.lbl, cr.lbl = cr.lbl,
                 prev = probs[1],
                 sens = probs[2],
                 spec = probs[4], fart = probs[5],
                 N = N,
                 hi = freqs$hi, mi = freqs$mi,
                 fa = freqs$fa, cr = freqs$cr,
                 scen.src = scen.src, scen.apa = scen.apa)

  ## Add class riskyr:
  class(object) <- "riskyr"

  return(object)

}
