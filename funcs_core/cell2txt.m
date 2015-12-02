function cell2txt(Cell2Trans,TxtName)
%     TxtName=input('Please Enter the Result TXT File Name (e.g. res) : ','s');

% codes below from
% http://stackoverflow.com/questions/8565617/print-a-cell-array-as-txt-in-matlab
% modified by hscheng
    fid = fopen([TxtName,'.txt'], 'a');
    for ii = 1:size(Cell2Trans,1)
        Temp=Cell2Trans{ii,1};
        fstr = '';
        
        % generate a print format
        for jj=1:size(Temp,2)
           
           switch class(Temp{jj})
               case 'char'
                   fstr = [fstr '%s'];
               otherwise
                   % Assume numeric
                   fstr = [fstr '%g'];
           end
           
           if jj < size(Temp,2)
               fstr = [fstr '\t'];
           else
               fstr = [fstr '\n'];
           end
        
        end
        
        C = Temp.';

        fprintf(fid, fstr,C{:});
    end

    fclose(fid);
end