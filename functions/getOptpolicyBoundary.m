function optbound = getOptpolicyBoundary(optpolicy,params)
% Get optimal policy boundaries from 3D optimal policy matrix
% Outputs the x,y,z values of the 3D boundaries
optbound = struct;

val = params.Rs';

fNames = {'dec','switch','accum'};
for fi = 1:length(fNames), optbound.(fNames{fi}) = cell(1,2); end
for y = 1:2
    for v = val
        val_i = val==v;
        thisslice = optpolicy.(sprintf('p%d',y))(:,:,val_i);
        for p = 1:4
            if any(thisslice(:)==p)
                thisslicemod = thisslice;
                thisslicemod(thisslicemod==p) = 99;
                thisslicemod(thisslicemod~=99) = 0;
                thisslicemod(thisslicemod==99) = 1;
                B = bwboundaries(thisslicemod);
                % remove edges
                thisBoundary = B{1};
                outerEdge_i = thisBoundary(:,1)==1 | thisBoundary(:,1)==length(params.ts) | thisBoundary(:,2)==1 | thisBoundary(:,2)==length(params.ts);
                thisBoundary(outerEdge_i,:) = [];
                % store the x,y coordinates and the r value
                if p==1 || p==2
                    optbound.dec{y} = cat(1,optbound.dec{y},[fliplr(thisBoundary),repmat(v,size(thisBoundary,1),1)]);
                elseif p==3
                    optbound.accum{y} = cat(1,optbound.accum{y},[fliplr(thisBoundary),repmat(v,size(thisBoundary,1),1)]);
                elseif p==4
                    optbound.switch{y} = cat(1,optbound.switch{y},[fliplr(thisBoundary),repmat(v,size(thisBoundary,1),1)]);
                end
            end
        end
    end
end

end

