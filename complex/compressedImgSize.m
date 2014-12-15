function sz = compressedImgSize(im, level)
import com.mathworks.mlwidgets.io.InterruptibleStreamCopier

if not(islogical(im))
  error('The image must be logical!!!!');
end

if isempty(im)
  sz = 0;
end

level = int32(level);
if (level<1) || (level>9)
  error('Level not correct %s!!!!!', any2str(level));
end


im = packLogicalImage(im);
im = im(:);

d=java.util.zip.Deflater(level);
f=java.io.ByteArrayOutputStream();
g=java.util.zip.DeflaterOutputStream(f,d);
g.write(im);
g.close();
d.end();
sz = f.size();
% Z=typecast(f.toByteArray,'uint8');
f.close;

% a=java.io.ByteArrayInputStream(Z);
% b=java.util.zip.InflaterInputStream(a);
% isc = InterruptibleStreamCopier.getInterruptibleStreamCopier;
% c = java.io.ByteArrayOutputStream;
% isc.copyStream(b,c);
% Q=typecast(c.toByteArray,'uint8');
% 
% Q=Q;
% 
% return
% 
