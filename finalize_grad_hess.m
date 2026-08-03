function [grad,hess]=finalize_grad_hess(raw_grad,raw_hess,want_hess,nx)

% FINALIZE_GRAD_HESS Squeeze raw (trailing-nx) grad/Hessian structs from
% BOUNDARY_RAW/RESOLVE_CUSP_TERMS into the public CDF_GRAD_BD/
% NORM_ERR_GRAD_BD return shapes.
%
% Internal helper, not meant to be called directly.
%
% Abhranil Das
% Center for Perceptual Systems, University of Texas at Austin
% Comments, questions, bugs to abhranil.das@utexas.edu

grad=struct();
if isfield(raw_grad,'q0')
    g0=raw_grad.q0; grad.q0=g0; if nx==1, grad.q0=g0(1); end
end
if isfield(raw_grad,'q1')
    Gq1=raw_grad.q1; grad.q1=Gq1; if nx==1, grad.q1=Gq1(:,1); end
end
if isfield(raw_grad,'q2')
    Gq2=raw_grad.q2; grad.q2=Gq2; if nx==1, grad.q2=Gq2(:,:,1); end
end

hess=[];
if want_hess
    hess=raw_hess;
    if nx==1                                     % drop the trailing singleton
        hess.q0q0=hess.q0q0(1); hess.q0q1=hess.q0q1(:,1);
        hess.q0q2=hess.q0q2(:,:,1); hess.q1q1=hess.q1q1(:,:,1);
        hess.q1q2=hess.q1q2(:,:,:,1); hess.q2q2=hess.q2q2(:,:,:,:,1);
    end
end

end
