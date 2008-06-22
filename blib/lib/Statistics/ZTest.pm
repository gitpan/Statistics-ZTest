package Statistics::ZTest;

use 5.008008;
use strict;
use warnings;
use Carp qw(croak);
use Class::OOorNO qw(coerce_array);
use vars qw($VERSION);
$VERSION = 0.01;

use Statistics::Distributions;

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sub new {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    my $class = shift;
    my $args = Class::OOorNO::coerce_array(@_);
	my $self = {};
	bless $self, $class;

    # Set default values:
    $self->{'tails'} = 2;
    $self->{'s_precision'} = 2;
    $self->{'p_precision'} = 0;
    $self->{'ccorr'} = 0;
    ##$self->{'distribution'} = 'norm';
    
    if (scalar keys %{$args}) {
        foreach (keys %{$args}) {
            $self->{$_} = $args->{$_};
        }
    }

	return $self;
}

#-----------------------------------------------
sub ztest {
#-----------------------------------------------
    my $self = shift;
    my $args = Class::OOorNO::coerce_array(@_);
    foreach (keys %{$args}) {
        $self->{$_} = $args->{$_};
    }

    my ($z, $p, $obs_dev, $std_dev) = ();

    $obs_dev = $self->{'observed'} - $self->{'expected'};
    $obs_dev = _continuity_correct($obs_dev) if $self->{'ccorr'} && $obs_dev;
    $std_dev = $self->{'std_dev'} ? 
        $self->{'std_dev'} : 
            $self->{'variance'} ?
        sqrt($self->{'variance'}) :
            return undef; #croak "Need a std_dev or variance value";
    $z = $obs_dev / $std_dev;
 
    $p = $self->p_value($z);
    $z = sprintf('%.' . $self->{'s_precision'} . 'f', $z) if $self->{'s_precision'};
    
    return wantarray ? ($z, $p, $obs_dev, $std_dev) : $z;
}

#-----------------------------------------------
sub p_value {
#-----------------------------------------------
    my ($self, $z) = @_;
    my $p = Statistics::Distributions::uprob(abs($z));
    $p *= 2 if ! $self->{'tails'} or $self->{'tails'} == 2 || ($self->{'tails'} != 1 && $self->{'tails'} != 2);
	$p = Statistics::Distributions::precision_string($p);
    $p = 1 if $p > 1;
    $p = 0 if $p < 0;
    $p = sprintf('%.' . $self->{'p_precision'} . 'f', $p) if $self->{'p_precision'};
    return $p;
}

#-----------------------------------------------
sub r_2_z {
#-----------------------------------------------
    my ($self, $r) = @_;
    croak 'Null or invalid r_value' if !defined $r or $r > 1 or $r < -1;
    my $z = .5 * ( log( (1 + $r) / (1 - $r) ) );
    return $z;
}

#-----------------------------------------------
sub z_2_r {
#-----------------------------------------------
   my ($self, $z) = @_;
   return defined $z ? (exp(2 * $z) - 1) / (exp(2 * $z) + 1) : undef;
}

#-----------------------------------------------
sub chi_2_z {
#-----------------------------------------------
    my ($self, $chi) = @_;
    return sqrt($chi);
}

#-----------------------------------------------
sub intercorrelation {
#-----------------------------------------------
    my ($self, $r12, $r13, $r23, $n) = @_;
    my $co_prod = ( $r12**2 + $r13**2 ) / 2;
    my $f_val = ( 1 - $r23 ) / ( 2 * (1 - $co_prod ) );
    my $h_val = ( 1 - $f_val * $co_prod ) / ( 1 - $co_prod );
    my $z_val = ( $self->r_2_z($r12) - $self->r_2_z($r13) ) * sqrt( ($n - 3) / ( 2 * (1 - $r23) * $h_val ) );
    return wantarray ? ($z_val, $self->p_value($z_val)) : $z_val;
}

#-----------------------------------------------
# Apply the continuity correction to the deviation, e.g. of (x - MCE) in numerator of z-score, if variance is calculated binomially (as hit/miss, not continuous):
sub _continuity_correct {
#-----------------------------------------------
   my $dev = shift;
   my $c = abs($dev) - .5;
   $c *= -1 if $dev < 0;
   return $c;
}

#-----------------------------------------------
sub state_significance {
#-----------------------------------------------
	my ($p, $dev) = @_;
	return 'The observed deviation was ' .
        (
           $p < .075 ? 
            ( 
               $p < .05 ? 
                    'statistically significant at the ' . 
                     (
                         $p < .01 ? '.01' : '.05'
                     ) . " level.\nThe observed value was " . 
                     (
                          $dev > 0 ? 'greater' : 'less' 
                      ) . 
                     ' than expected by chance' :
                  'of marginal statistical significance'
               ) :
               'not statistically significant'
           ) .
    ".\n";
}


# --------------------
# Function aliases
# --------------------
*test = \&ztest;
*zvalue = \&ztest;
*meng_intercorrelation = \&intercorrelation;

1;
__END__

=head1 NAME

Statistics::ZTest - Basic ztest and normal probability reporting

=head1 VERSION

This is documentation for Version 0.01 of Statistics::ZTest (2006.11.22).

=head1 SYNOPSIS

 use Statistics::ZTest;
 my $dev = Statistics::ZTest->new(
    ccorr    => 1,
    tails    => 2,
    s_precision => 5,
    p_precision => 5,
 );

 my ($z, $pz, $observed_deviation, $standard_deviation) = 
    $dev->ztest(
        observed => $stat_obs,
        expected => $stat_exp,
        variance => $variance,
 );

=head1 DESCRIPTION

Calculates a z-statistic: the ratio of an observed deviation to a standard deviation. Purpose is simply to support L<Statistics::Sequences|Statistics::Sequences>, but with some standalone utilities.

=head1 METHODS

=head2 new

 $dev = Statistics::ZTest->new();

Returns a Statistics::ZTest object. Accepts setting of any of the L<OPTIONS>.

=head2 ztest

You supply the observed and expected values of your statistic, and the variance (observed or expected).

Additionally, you may specify a boolean for performing the continuity-correction to the observed deviation, and then a value of either 1 or 2 to specify the tails relevant to determining the probability of obtaining the calculated z.

Returns an array consisting of the z-statistic, its probability, the observed deviation (the difference between the observed and expected values of your statistic), and the standard deviation (the square-root of the variance supplied).

All return a z_value, or a z_value and p_value, where relevant, and called in array context.

=head2 p_value

 $p = $dev->p_value($z)

Send a z-value, get its associated p-value, 2-tailed by default.

=head2 r_2_z

 $z = $dev->r_2_z($r)

Performs the Fisher r-to-z transformation. Send a correlation coefficient - get back a z-value.

=head2 z_2_r

 $r = $dev->z_2_r($z)

Send a z-value - get back a correlation coefficient.

=head2 chi_2_z

 $z = $dev->chi_2_z($chi) 

Send a chi-value, get back a z-value (the square-root of the thing ...).

=head2 dep_intercorrelation

 ($z, $p) = $dev->dep_intercorrelation($r_xz, $r_yz, $r_xy, $n)

For comparing dependent correlations with a common variable. You have the r-values for the correlation between (1) Variable X and Variable Z, and Variable Y and Variable Z, and want to test, for example, if Z correlates more highly with X than with Y. These two correlations form the first two arguments (after the class object). You also need to calculate the correlation of X and Y themselves, and send it as the third argument. Lastly, also send the sample size. Uses the Meng-Rosenthal-Rubin method.

=head1 OPTIONS

The following can be set in the call to L<new> or L<test>.

=head2 observed

The observed value of the test statistic.

=head2 expected

The expected value of the test statistic.

=head2 variance

The variance of the test statistic (whether expected or observed).

=head2 ccorr

Apply the continuity correction. Default = 0.

=head2 s_precision

Precision of the z_value. Default = 2.

=head2 p_precision

Precision of the associated p_value. Default = 0.

=head2 tails

Tails from which to assess the association p_value (1 or 2). Default = 2.

=head1 SEE ALSO

L<Statistics::Distributions|Statistics::Distributions> : the C<uprob> and C<precision_string> methods are here used for calculating and reporting probability.

=head1 TO DO/BUGS

Other distributions.

=head1 AUTHOR

Roderick Garton, E<lt>rgarton@utas_DOT_edu_DOT_auE<gt>

=head1 COPYRIGHT/LICENSE/DISCLAIMER

Copyright (C) 2007 Roderick Garton 

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself, either Perl version 5.8.8 or, at your option, any later version of Perl 5 you may have available. 

To the maximum extent permitted by applicable law, the author of this module disclaims all warranties, either express or implied, including but not limited to implied warranties of merchantability and fitness for a particular purpose, with regard to the software and the accompanying documentation.

=cut


1;
