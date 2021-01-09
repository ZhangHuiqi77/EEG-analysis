function [] = speakFunc(name,content)
    if ~exist('content','var')
        content = 'hello MATLAB';
    end
    fprintf(1,'%s speak: %s\n',name,content);
end