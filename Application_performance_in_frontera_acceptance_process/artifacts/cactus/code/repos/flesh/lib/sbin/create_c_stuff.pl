#! /usr/bin/perl
#/*@@
#  @file      create_c_stuff.pl
#  @date      Mon Jan 11 10:53:22 1999
#  @author    Tom Goodale
#  @desc
#
#  @enddesc
#  @version $Id$
#@@*/


#/*@@
#  @routine    CreateParameterBindingFile
#  @date       Wed Jan 20 15:20:23 1999
#  @author     Tom Goodale
#  @desc
#  Creates the bindings used to link the thorn parameters with the flesh.
#  @enddesc
#@@*/

sub CreateParameterBindingFile
{
  my($include_headers, $structure, $rhparameters, $rhparameter_db) = @_;
  my($line,@data);
  my(%parameters);
  my($type);

  # Header Data
  if ($include_headers)
  {
    push(@data, '#include "cctk_Config.h"');
    push(@data, '#include "CParameterStructNames.h"');
    push(@data, '');
  }

  # Create the structure
  push(@data, 'struct');
  push(@data, '{');

  foreach $parameter (&order_params($rhparameters,$rhparameter_db))
  {
    my $type = $rhparameter_db->{"\U$rhparameters->{$parameter} $parameter\E type"};
    my $type_string = &get_c_type_string($type,$parameter);

    my $realname = $rhparameter_db->{"\U$rhparameters->{$parameter} $parameter\E realname"};

    my $array_size = $rhparameter_db->{"\U$rhparameters->{$parameter} $parameter\E array_size"};

    my $suffix = '';

    if($array_size)
    {
      $suffix = "[$array_size]";
    }

    push(@data, "  $type_string $realname$suffix;");
  }

  # Some compilers don't like an empty structure.
  if((keys %$rhparameters) == 0)
  {
    push(@data, '  int dummy_parameter;');
  }

  push(@data, "} $structure;");
  push(@data, "\n");   # workaround for perl 5.004_04 to add a trailing newline

  return join ("\n", @data);
}


#/*@@
#  @routine    get_c_type_string
#  @date       Mon Jan 11 15:33:50 1999
#  @author     Tom Goodale
#  @desc
#  Returns the correct type string for a parameter.
#  @enddesc
#@@*/

sub get_c_type_string
{
  my($type,$parameter) = @_;
  my $type_string;


  if($type eq 'KEYWORD' ||
     $type eq 'STRING'  ||
     $type eq 'SENTENCE')
  {
    $type_string = 'const char *';
  }
  elsif($type eq 'BOOLEAN')
  {
    $type_string = 'CCTK_INT';
  }
  elsif($type eq 'INT')
  {
    $type_string = 'CCTK_INT';
  }
  elsif($type eq 'INT1')
  {
    $type_string = 'CCTK_INT1';
  }
  elsif($type eq 'INT2')
  {
    $type_string = 'CCTK_INT2';
  }
  elsif($type eq 'INT4')
  {
    $type_string = 'CCTK_INT4';
  }
  elsif($type eq 'INT8')
  {
    $type_string = 'CCTK_INT8';
  }
  elsif($type eq 'INT16')
  {
    $type_string = 'CCTK_INT16';
  }
  elsif($type eq 'REAL')
  {
    $type_string = 'CCTK_REAL';
  }
  elsif($type eq 'REAL4')
  {
    $type_string = 'CCTK_REAL4';
  }
  elsif($type eq 'REAL8')
  {
    $type_string = 'CCTK_REAL8';
  }
  elsif($type eq 'REAL16')
  {
    $type_string = 'CCTK_REAL16';
  }
  else
  {
    &CST_error(0,"Unknown parameter type '$type' for parameter '$parameter'",'',__LINE__,__FILE__);
  }

  return $type_string;
}


#/*@@
#  @routine    GetThornParameterList
#  @date       Wed Jan 20 15:29:40 1999
#  @author     Tom Goodale
#  @desc
#  Gets a list of all parameters in a particular block in a thorn.
#  Returns a hash table.
#  @enddesc
#@@*/

sub GetThornParameterList
{
  my($thorn, $block, $rhparameter_db) = @_;
  my(%parameter_list);

  $params = $rhparameter_db->{"\U$thorn $block\E variables"};

  foreach $parameter (split(' ', $params))
  {
    if($parameter =~ m:[^ ]:)
    {
      $parameter_list{$parameter} = $thorn;
    }
  }

  return %parameter_list;
}


sub CreateCStructureParameterHeader
{
  my($prefix, $structure, $rhparameters, $rhparameter_db) = @_;
  my($line,@data);
  my(%parameters);
  my($type, $type_string);
  my(@definition, @definition2);

  # determine thorn name from structure name
  # TODO: pass thorn name explicitly
  my $thorn = $structure;
  $thorn =~ s{^[^_]*_}{};       # remove PRIVATE_ or RESTRICTED_ prefix
  $thorn =~ s{_STRUCT$}{};

  # Create the structure
  push(@data, '#ifdef __cplusplus');
  push(@data, 'extern "C"');
  push(@data, '{');
  push(@data, '#endif');
  push(@data, '');
  push(@data, 'extern struct');
  push(@data, '{');

  foreach $parameter (&order_params($rhparameters, $rhparameter_db))
  {
    my $type = $rhparameter_db->{"\U$rhparameters->{$parameter} $parameter\E type"};
    my $type_string = &get_c_type_string($type);

    my $array_size = $rhparameter_db->{"\U$rhparameters->{$parameter} $parameter\E array_size"};

    my $suffix = '';
    my $varprefix = '';

    if($array_size)
    {
      $varprefix = ' const *';
      $suffix = "[$array_size]";
    }

    my $realname = $rhparameter_db->{"\U$rhparameters->{$parameter} $parameter\E realname"};

    push(@data, "  $type_string $realname$suffix;");
    #push(@definition, "  CCTK_DECLARE_INIT ($type_string$varprefix const, $parameter, $structure.$realname); \\");
    push(@definition, "  CCTK_DECLARE_INIT ($type_string$varprefix const, $parameter, CCTK_PARAMETER__${thorn}__$realname); \\");
    push(@definition2,
         "#ifndef CCTK_PARAMETER__${thorn}__$realname\n" .
         "#  define CCTK_PARAMETER__${thorn}__$realname $structure.$realname\n" .
         "#endif");
  }

  # Some compilers don't like an empty structure.
  if((keys %$rhparameters) == 0)
  {
    push(@data, '  int dummy_parameter;');
  }

  push(@data, "} $structure;");
  push(@data, '');

  push(@data, '#ifdef __cplusplus');
  push(@data, '}');
  push(@data, '#endif');
  push(@data, '');

  push(@data, "#define DECLARE_${structure}_PARAMS \\");
  push(@data, @definition);
  push(@data, "\n");
  push(@data, @definition2);
  push(@data, "\n");

  return join ("\n", @data);
}


sub order_params
{
  my($rhparameters, $rhparameter_db) = @_;
  my(@float_params) = ();;
  my(@int_params)   = ();
  my(@string_params)= ();

  foreach $parameter (sort(keys %$rhparameters))
  {
    $type = $rhparameter_db->{"\U$rhparameters->{$parameter} $parameter\E type"};

    if($type eq 'KEYWORD' ||
       $type eq 'STRING'  ||
       $type eq 'SENTENCE')
    {
      push(@string_params, $parameter);
    }
    elsif($type eq 'BOOLEAN' ||
          $type eq 'INT')
    {
      push(@int_params, $parameter);
    }
    elsif($type eq 'REAL')
    {
      push(@float_params, $parameter);
    }
    else
    {
      $message = "Unknown parameter type '$type'";
      &CST_error(0,$message,__LINE__,__FILE__);
    }
  }

  return (@float_params, @string_params, @int_params);
}


#sub create_parameter_code
#{
#  my($structure, $implementation,$parameter, $rhparameter_db) = @_;
#  my($type, $type_string);
#  my($line, @lines);
#  my($default);
#  my($temp_default);
#
#  $default = $rhparameter_db->{"\U$implementation $parameter\E default"};
#  $type = $rhparameter_db->{"\U$implementation $parameter\E type"};
#
#  $type_string = &get_c_type_string($type);
#
#  if($type_string eq 'char *')
#  {
#    $line = "  $structure.$parameter = malloc(" . (length($default)-1) . '*sizeof(char));';
#    push(@lines, $line);
#
#    push(@lines, "  if ($structure.$parameter)");
#    push(@lines, "    strcpy($structure.$parameter, $default);");
#  }
#  elsif($type eq "BOOLEAN")
#  {
#    # Logicals need to be done specially.
#
#    # Strip out any quote marks, and spaces at start and end.
#    $temp_default = $default;
#    $temp_default =~ s:\"::g;
#    $temp_default =~ s:\s*$:: ;
#    $temp_default =~ s:^\s*:: ;
#
#    push(@lines, "  CCTK_SetLogical(\&($structure.$parameter),\"$temp_default\");");
#  }
#  else
#  {
#    push(@lines, "  $structure.$parameter = $default;");
#  }
#
#  $line = "CCTKi_ParameterCreate($parameter, $implementation,
#                    \"foobar\",\"" . $rhparameter_db->{"\U$implementation $parameter\E type"}."\"
#                    const char *scope,
#                    int        steerable,
#                    const char *description,
#                    const char *defval,
#                    void       *data)";
#
#  return @lines;
#}

1;
