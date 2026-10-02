function data = read_binary_img(filename)
    % 检查输入参数
    if nargin < 1
        error('必须指定要读取的文件名');
    end
    
    % 打开二进制文件（以二进制读取模式）
    fid = fopen(filename, 'rb');
    if fid == -1
        error('无法打开文件 "%s"。请检查文件路径和权限。', filename);
    end
    
    % 读取整个文件内容为single类型（float32）的列向量
    data = fread(fid, inf, 'single');
    
    % 关闭文件
    fclose(fid);
    
    % 确保返回的数据是single类型
    data = single(data);
end