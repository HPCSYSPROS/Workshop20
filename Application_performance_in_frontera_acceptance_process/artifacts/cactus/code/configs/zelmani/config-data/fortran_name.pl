
sub fortran_name
{
  my($old_name) = @_;
  my($new_name);

  $new_name = "\L$old_name\E";

  if($new_name =~ m:_: )
  {
    $new_name = $new_name."_";
  }
  else
  {
    $new_name = $new_name."_";
  }

  return $new_name;
}


sub fortran_common_name
{
  my ($old_name) = @_;
  my ($new_name);

  $new_name = "\L$old_name\E";
  if($new_name =~ m:_: )
  {
      $new_name = $new_name."_";
  }
  else
  {
      $new_name = $new_name."_";
  }

  return "".$new_name;
}

1;
