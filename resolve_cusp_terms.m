function [resolved_grad,resolved_hess]=resolve_cusp_terms(raw_grad,raw_hess,at_cusp,inversion_grad_hess,cusp_scale,nx)

% RESOLVE_CUSP_TERMS Override only the at_cusp points of raw_grad/raw_hess
% (in the same raw, trailing-nx shape as BOUNDARY_RAW returns), leaving
% every other point -- and every other problem that never hits this
% measure-adjacent exact coincidence -- completely untouched.
%
% inversion_grad_hess(s_eff,want_hess) is the probe: the (always valid,
% cusp-independent) s~=0 inversion route at a shrinking sequence of
% artificial normal terms, in the same raw shape as raw_grad/raw_hess.
% CDF_GRAD_BD passes a single class's own probe; NORM_ERR_GRAD_BD
% passes the *combined* p1*F1-p0*F0 probe instead, so that a real
% cancellation between the two classes' divergences shows up as convergence
% here rather than being lost by resolving each class's infinity
% independently beforehand.
%
% Internal helper, not meant to be called directly.
%
% Abhranil Das
% Center for Perceptual Systems, University of Texas at Austin
% Comments, questions, bugs to abhranil.das@utexas.edu

% Shrinking artificial normal-term probes used to test convergence at the
% cusp. Chosen from direct testing on the term that first exposed this
% (problem 6 of the wider benchmark): a genuinely divergent term shows a
% clean, unmistakable, non-shrinking decade-over-decade change throughout
% this whole range (and well beyond it); a genuinely finite term converges
% within the first one or two decades of it. Below roughly 1e-8*scale the
% inversion route itself becomes ill-conditioned (checked directly)
% regardless of convergence, so this stays well clear of that.
CUSP_PROBE_FACTORS=[1e-2 1e-3 1e-4];

want_hess=~isempty(raw_hess);
np_=numel(CUSP_PROBE_FACTORS);
probes_g=cell(1,np_); probes_h=cell(1,np_);
for pp=1:np_
    [gp,hp]=inversion_grad_hess(CUSP_PROBE_FACTORS(pp)*cusp_scale,want_hess);
    probes_g{pp}=gp;
    if want_hess, probes_h{pp}=hp; end
end

resolved_grad=raw_grad;
gfields=fieldnames(raw_grad);
for fi=1:numel(gfields)
    ky=gfields{fi};
    merged=cusp_probe_merge({probes_g{1}.(ky),probes_g{2}.(ky),probes_g{3}.(ky)});
    resolved_grad.(ky)=cusp_select(merged,raw_grad.(ky),at_cusp,nx);
end

if ~want_hess, resolved_hess=[]; return; end
resolved_hess=raw_hess;
hfields=fieldnames(raw_hess);
for fi=1:numel(hfields)
    ky=hfields{fi};
    merged=cusp_probe_merge({probes_h{1}.(ky),probes_h{2}.(ky),probes_h{3}.(ky)});
    resolved_hess.(ky)=cusp_select(merged,raw_hess.(ky),at_cusp,nx);
end

end

function y=cusp_select(merged,original,at_cusp,nx)
% elementwise select merged where at_cusp broadcasts true against the
% trailing nx axis, else original. nx==1 is special-cased: MATLAB silently
% drops a trailing singleton dimension (size(reshape(X,d,d,1)) is [d,d],
% not [d,d,1]), so inferring the nx axis from merged's shape would
% misidentify the last "d" axis as nx and produce a wrong, undersized mask
% -- at_cusp is then just a single scalar, so no per-element broadcasting
% is needed at all.
if nx==1
    if at_cusp, y=merged; else, y=original; end
    return
end
sz=size(merged);
shp=ones(1,numel(sz)); shp(end)=nx;
mask=repmat(reshape(at_cusp,shp),[sz(1:end-1),1]);
y=original;
y(mask)=merged(mask);
end

function out=cusp_probe_merge(vals)
% vals: 1-by-3 cell of same-shape arrays evaluated at a shrinking sequence
% of tiny artificial normal-term probes (largest first). Elementwise: if
% converged (latest decade-over-decade change well under the one before),
% trust the smallest-probe value; otherwise the term is genuinely unbounded
% and the correctly-signed infinity is returned instead of an arbitrary
% finite number.
v0=vals{1}; v1=vals{2}; v2=vals{3};
d1=v1-v0; d2=v2-v1;
converging=(d1==0)|(abs(d2)<0.5*abs(d1));
out=v2;
div=~converging;
out(div & (d2>0))=inf;
out(div & ~(d2>0))=-inf;
end
