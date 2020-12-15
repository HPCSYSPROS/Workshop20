#! /usr/bin/perl -s
#/*@@
#  @file      Logger
#  @date      Thu April 19 14:25:13 2004
#  @author    Yaakoub El-Khamra
#  @desc
#             Logs configuration scripts
#  @enddesc
#  @version   $Header$
#@@*/


require "$sbin_dir/CSTUtils.pl";

#/*@@
#  @routine   CreateLogFile
#  @date
#  @author    Yaakoub El-Khamra
#  @desc
#  Writes logs of configuration scripts
#  @enddesc
#@@*/
sub CreateLogFile
{
  my($config_dir, $cfg, $thorns)=@_;
  my($data) = '';

  foreach $thorn (sort keys %thorns)
  {
    if ($cfg->{"\U$thorn\E PROVIDES"})
    {
      foreach $provides (split (' ', $cfg->{"\U$thorn\E PROVIDES"}))
      {
        $data .= "Configuration Script log for thorn: $thorn:\n";
        $data .= $cfg->{"\U$thorn $provides\E MESSAGE"};
        $data .= $cfg->{"\U$thorn $provides\E ERROR"};
      }
    }
  }

  &WriteFile("$config_dir/ConfigScripts.log", \$data);

}


1;
