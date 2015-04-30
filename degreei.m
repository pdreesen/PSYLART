function degi = degreei(term, vari),
% degreei(term, vari) 
% computes degree of variable 'vari' in term 'term'

if length(term) > 1,
    warning( 'several terms: considering only first term' );
    term = term(1);
end

degi = double( coeffs( diff( term, vari ) ) ); 

if isempty(degi), degi = 0; end


end