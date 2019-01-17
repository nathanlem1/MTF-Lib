function keys= unique_faster(keys)

keys= sort(keys(:));
difference = diff([keys;NaN]);
keys = keys(difference~=0)';

end
